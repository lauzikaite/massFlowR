#' @title Group peaks of a single LC-MS spectrum using EIC correlation
#' @description Function builds peak-groups of co-eluting peaks that are correlated to each other in a single LC-MS spectrum.
#'
#' @param pks A \code{DataFrame} class peak table, created by \emph{pickPEAKS} function.
#' @param eic A \code{list} containing extracted ion chromatograms for each peak in the \code{pks} table.
#' @param match A \code{numeric} defining the number of scans for co-eluting peaks extraction.
#' @param thr A \code{numeric} defining correlation coefficient threshold, above which peak pairs will be considered as correlated.
#' @param plot A \code{logical}. For \code{plot = TRUE}, a network graph for each peak in the table will be saved in the out_dir directory.
#' @param out_dir A \code{character} specifying directory where output data will be saved.
#' @param fname A \code{character} specifying LC-MS filename.
#' @param pearson A \code{logical} whether Pearson Correlation should be used. For \code{pearson = FALSE}, Spearman correlation method will be used.
#' @param clean A \code{logical} whether one-peak peak-groups should be removed (default is TRUE).
#'
#' @return Function returns an updated peak table. Column \code{"peakgr"} contains peak-group ID for each peak.
#' Peaks not correlated to any of its co-eluting peaks are removed from the table.
#' Function also writes updated peak table in the specified directory.
#' 
#' @export
#'

groupPEAKSspec <- function(pks, eic, out_dir, fname, pearson = TRUE, match = 1, thr = 0.95, plot = FALSE, clean = TRUE, return = FALSE) {

  message("Building peak-groups for file: ", fname, "...")

  ## make a table to store peak-groups during grouping
  peakgroups <- pks %>%
    mutate(peakgr = NA) %>%
    select(peakid, mz, scpos, into, peakgr) %>%
    data.frame()

  peakids <- peakgroups %>%
    pull(peakid)

  ####--- (A) select POI (peak of interest) from the table
  for (p in peakids) {
    if (!is.na(peakgroups[p, "peakgr"])) next ## if POI is already assigned to a group, skip to next

    ## find co-eluting peaks not asigned to a peakgronent yet
    if (match == 0) {
      ## exact scpos match
      pks_co <- peakgroups %>%
        filter(is.na(peakgr) & scpos == scpos[peakid == p])
    } else {
      ## scpos +- match
      mat <- c(peakgroups[p, "scpos"] - match, peakgroups[p, "scpos"] + match)
      pks_co <- peakgroups %>%
        filter(is.na(peakgr)) %>% filter(between(scpos, mat[1], mat[2]))
    }
    ## store indeces of co-eluting peaks
    co_ind <- pks_co %>%
      pull(peakid)
    if (length(co_ind) == 1 | !p %in% co_ind) next ## if POI doesn't have co-eluting peaks, skip to next

    ####--- (B) build correlation matrix for all co-eluting peaks
    poi_co_cor <- buildCOR(co_ind = co_ind, eic = eic, pearson = pearson)
    if (all(poi_co_cor$cor < thr)) next ## if there isn't a single pair with cor above threshold, skip to next

    ####---- (C) build interaction network between co-eluting peaks
    ## returns updated table with assigned peak-group IDs
    peakgroups <- buildNETWORK(poi_co_cor = poi_co_cor,
                             peakgroups = peakgroups,
                             co_ind = co_ind,
                             thr = thr,
                             pks = pks,
                             p = p,
                             plot = plot,
                             out_dir = out_dir,
                             fname = fname)

  }

  if (clean == TRUE) {
    ## removing peaks not assigned to any peak-group
    peakgroups <- peakgroups %>% filter(!is.na(peakgr))
  } else {
    peakgroups[which(is.na(peakgroups$peakgr)), "peakgr"] <- seq(
        from = (max(na.omit(peakgroups$peakgr)) + 1),
        length.out = length(which(is.na(peakgroups$peakgr))),
        by = 1)
  }

  peakgroups <- merge(pks, peakgroups[, c("peakid", "peakgr")],
                    by = c("peakid"), all = F)
  write.csv(peakgroups, file = paste0(out_dir, "/", fname, "_peakgrs.csv"), quote = F, row.names = F)

  if (return == TRUE) {
    return(peakgroups)
  } else {
    message(max(na.omit(peakgroups$peakgr)), " peak-groups built.")
  }

}


####---- helper functions, not to export ----
## Build correlation matrix using full EIC list and stating which peaks from that list to correlate
buildCOR <- function(co_ind, eic, pearson) {
  poi_co_cor <- expand.grid(x = co_ind, y = co_ind) %>%
    group_by(x, y) %>%
    mutate(cor = getCOR(x = x, y = y, eic = eic, pearson = pearson)) %>%
    filter(x != y)
}

## Correlate two peaks using correspondig EICs from the list of EIC
getCOR <- function(x, y, eic, pearson) {
  rx <- MSnbase::rtime(eic[[x]])
  ry <- MSnbase::rtime(eic[[y]])
  common_scan <- base::intersect(rx, ry)
  if (length(common_scan) > 3) {
    ix <- as.numeric(MSnbase::intensity(eic[[x]])[which(rx %in% common_scan)])
    iy <- as.numeric(MSnbase::intensity(eic[[y]])[which(ry %in% common_scan)])
    cc <- ifelse(pearson == FALSE, cor(ix, iy, method = "spearman"), cor(ix, iy, method = "pearson"))
  } else {
    cc <- 0
  }
  return(cc)
}

## build network using correlation coefficients
buildNETWORK <- function(poi_co_cor, peakgroups, co_ind, match, thr, pks, p, plot, out_dir, fname, return = TRUE, colors = NULL) {

  cormat <- poi_co_cor %>%
    tidyr::spread(y, cor) %>%
    data.frame
  rownames(cormat) <- cormat$x
  cormat$x <- NULL
  cormat <- as.matrix(cormat)

  ## build network with all feature pairs, ignore self-self correlation
  g <- igraph::graph.adjacency(cormat, weighted = TRUE, mode = "lower", diag = FALSE)

  ## asign names to vertices(nodes)
  igraph::V(g)$name <- peakgroups %>%
    filter(peakid %in% co_ind) %>%
    pull(mz) %>%
    round(digits = 3)

  ## delete edges between feature pairs with cor < thr
  g <- igraph::delete.edges(g, igraph::E(g)[[which(igraph::E(g)$weight < thr)]])

  ## detect community structure in the network where non-correlating pairs are omitted
  coms <- igraph::cluster_label_prop(g)

  ## color all nodes according to their status: 1 - peak-of-interest (POI), 2 - belongs to the same community as POI, 3 - unrelated
  mem <- igraph::membership(coms)
  poi <- mem[which(names(mem) == round(pks[p, "mz"], digits = 3))]
  mem <- ifelse(mem == poi, 2, 3)
  mem[which(names(mem) == round(pks[p, "mz"], digits = 3))] <- 1

  ## delete edges between different communities
  gg <- igraph::delete.edges(g, igraph::E(g)[igraph::crossing(coms, g)])

  if (plot == TRUE) plotNETWORK(gg = gg, mem = mem, out_dir = out_dir, match = match, fname = fname, p = p, co_ind = co_ind, colors = colors)

  ## if there is atleast one peak in the same community as the POI (i.e. in the mem == 2)
  members <- which(mem == 2 | mem == 1)

  if (length(members) >= 2) {
    ## which of the co-eluting peaks are part of this network community
    peak_ids <- co_ind[members]

    ## update table with assigned peakgroup IDs
    ## start with 1 if this is the first peakgroup to be assigned
    id <- ifelse(all(is.na(peakgroups$peakgr)), 1, max(peakgroups$peakgr, na.rm = T) + 1)
    peakgroups[peak_ids, "peakgr"] <- rep(id, length(peak_ids))
  }
  if (return == TRUE) { return(peakgroups) } else { message("Network built was succesfull")}
}

plotNETWORK <- function(gg, mem, out_dir, fname, match, p, co_ind, colors) {

  grid::grid.newpage()
  grDevices::pdf(width = 8, height = 8,
      file = paste0(out_dir, "/", fname, "_peak", p, "_CoNetwork.pdf"))

  ## Peak number of the peak of interest
  pmain <- paste0("Main peak: ", round(pks[p, "mz"], digits = 4),
                  ". Coeluting peaks: ",
                  length(co_ind))
  psub <- paste0("Features in POI peakgroup: ",
                 length(which(mem == 1 | mem == 2)))

  ## default color selection
  if (is.null(colors)) {
    colors <- c("#ABDDA4", "#E6F598", "grey")
    vertex.color = colors[mem]
  } else{
    vertex.color = colors[names(mem)]
  }

  ## make coordinates for vertices based on correlation: scale cor coef to maximise distance
  ## using Fruchterman-Reingold layout algorithm to prevent overlap
  coords <- igraph::layout_(gg, igraph::with_fr(weights = scaleEDGES(igraph::E(gg)$weight, from = 0.01, to = 10)))

  igraph::plot.igraph(gg,
                      main = pmain,
                      sub = psub,
                      layout = coords,
                      vertex.color = vertex.color,
                      vertex.size = 22,
                      vertex.label.color = "black",
                      edge.color = "royalblue4",
                      # edge.label = round(igraph::E(gg)$weight, digits = 2), ## labeling edges with cor makes it messy
                      edge.label.color = "royalblue4")
  grDevices::dev.off()
}

## function to scale correlation coefficients to make a more clear Fruchterman-Reingold graph
scaleEDGES <- function(x, from, to) {
  (x - min(x)) / max(x - min(x)) * (to - from) + from
  }


