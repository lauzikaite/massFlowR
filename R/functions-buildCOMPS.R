#' @title Component building based on correlation between co-eluting peaks.
#'
#' @param pks \code{DataFrame} object, provided by \emph{pickPEAKS} function.
#' @param eic \code{list} containing extracted ion chromatograms for each peak in the \code{pks} table.
#' @param match \code{numeric} defining the number of scans for co-eluting peaks extraction.
#' @param thr \code{numeric} defining correlation coefficient threshold, above which peak pairs will be considered as correlated.
#' @param plot \code{logical}. For \code{plot = TRUE}, a network graph for each peak in the table will be saved as a png in the out_dir directory.
#' @param out_dir \code{character} object specifying directory where output data will be saved.
#' @param fname \code{character} object specifying datafile name.
#' @param pearson \code{logical} whether Pearson Correlation should be used. For \code{pearson = FALSE}, Spearman correlation method will be used.
#' @param clean \code{logical} whether one-peak components should be removed (default is TRUE).
#'
#' @return Function returns a \code{DataFrame} object with a column \code{"comp"}. \code{"comp"} contains component ID for each peak, which was correlated with any of its co-eluting peak(s).
#' @export
#'
#' @examples

buildCOMPS <- function(pks, eic, out_dir, fname, pearson = TRUE, match = 1, thr = 0.95, plot = FALSE, clean = TRUE) {

  if (missing(fname)) { stop("'fname' must be specified!") }

  message("Building components for file:", fname, "...")

  ## duplicate table for storing built componenents
  pkscomps <- pks %>%
    mutate(comp = NA) %>%
    select(pid = pid, mz, scpos, into, comp) %>%
    data.frame

  pids <- pkscomps %>%
    pull(pid)

  # ## set progress bar
  # pb <- txtProgressBar(min = 0, max = nrow(pks), style = 3)

  ####--- (A) select POI (peak of interest) from the table ----
  for (p in pids) {
    if (is.na(pkscomps[p, "comp"])) {

      # ## update progress bar
      # setTxtProgressBar(pb = pb, value = p)

      ####---- (A1) find co-eluting peaks not asigned to a component yet ----
      if (match == 0) {
        ## exact scpos match
        pks_co <- pkscomps %>%
          filter(is.na(comp) & scpos == scpos[pid == p])

      } else {
        ## scpos +- match
        mat <- c(pkscomps[p, "scpos"] - match, pkscomps[p, "scpos"] + match)
        pks_co <- pkscomps %>%
          filter(is.na(comp)) %>% filter(between(scpos, mat[1], mat[2]))
      }

      ## store indeces of co-eluting peaks
      co_ind <- pks_co %>%
        pull(pid)

      ####--- if POI has any co-eluting peaks ----
      if (length(co_ind) > 1 & p %in% co_ind) {

        ####--- build correlation matrix for all co-eluting peaks ----

        poi_co_cor <- buildCOR(co_ind = co_ind, eic = eic, pearson = pearson)

        ####---- build interaction network between co-eluting peaks ----

        ## if there's at least one pair with cor above threshold
        if (any(poi_co_cor$cor > thr)) {

          pkscomps <- buildNETWORK(poi_co_cor = poi_co_cor,
                                   pkscomps = pkscomps,
                                   co_ind = co_ind,
                                   thr = thr,
                                   pks = pks,
                                   p = p,
                                   plot = plot,
                                   out_dir = out_dir,
                                   fname = fname)

          } else next
        } else next
      } else next
  }

  message(max(na.omit(pkscomps$comp)), " components built.")

  if (clean == TRUE) {
    message("'clean' set to TRUE. Removing un-grouped peaks ...")
    message(length(which(is.na(pkscomps$comp))), " peaks removed.")
    pkscomps <- pkscomps %>% filter(!is.na(comp))
  } else {
    message("'clean' set to FALSE. Returning un-grouped peaks as 1-peak components ...")
    message(length(which(is.na(pkscomps$comp))), " 1-peak components built.")
    pkscomps[which(is.na(pkscomps$comp)), "comp"] <- seq(
        from = (max(na.omit(pkscomps$comp)) + 1),
        length.out = length(which(is.na(pkscomps$comp))),
        by = 1)
  }

  pkscomps <- merge(pks, pkscomps[, c("pid", "comp")],
                    by = c("pid"), all = F)
  write.table(pkscomps, file = paste0(out_dir, fname, "_pks-comps.txt"), col.names = T, quote = F, sep = "\t", row.names = F)
  return(pkscomps)
}


####---- helper functions, not to export ----
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

buildNETWORK <- function(poi_co_cor, pkscomps, co_ind, match, thr, pks, p, plot, out_dir, fname, return = TRUE) {

  cormat <- poi_co_cor %>%
    tidyr::spread(y, cor) %>%
    data.frame

  rownames(cormat) <- cormat$x
  cormat$x <- NULL
  cormat <- as.matrix(cormat)

  ## build network with all feature pairs, ignore self-self correlation
  g <- igraph::graph.adjacency(cormat, weighted = TRUE, mode = "lower", diag = FALSE)

  ## asign names to vertices(nodes)
  igraph::V(g)$name <- pkscomps %>%
    filter(pid %in% co_ind) %>%
    pull(mz)

  ## delete edges between feature pairs with cor < thr
  g <- igraph::delete.edges(g, igraph::E(g)[[which(igraph::E(g)$weight < thr)]])

  ## detect community structure in the network where non-correlating pairs are omitted
  coms <- igraph::cluster_label_prop(g)

  ## color all nodes according to their status: 1 - POI, 2 - belongs to the same group as POI, 3 - unrelated
  mem <- igraph::membership(coms)
  poi <- mem[which(names(mem) == pks[p, "mz"])]
  mem <- ifelse(mem == poi, 2, 3)
  mem[which(names(mem) == pks[p, "mz"])] <- 1

  ## delete edges between different clusters
  gg <- igraph::delete.edges(g, igraph::E(g)[igraph::crossing(coms, g)])
  igraph::V(gg)$name <- pkscomps %>%
    filter(pid %in% co_ind) %>%
    pull(mz) %>%
    round(digits = 3)

  if (plot == TRUE) {
    plotNETWORK(gg = gg, mem = mem, out_dir = out_dir, match = match, fname = fname, p = p, co_ind = co_ind)
  }

  ## if there is at least once community with more than one member, asign community's compounds to a component
  members <- which(mem == 2 | mem == 1)

  if (length(members) >= 2) {
    component <- co_ind[members]

    ## update table with assigned component IDs
    comp <- ifelse(all(is.na(pkscomps$comp)), 1, max(pkscomps$comp, na.rm = T) + 1)
    pkscomps[component, "comp"] <- rep(comp, length(component))
  }
  if (return == TRUE) { return(pkscomps) } else { message("Network built was succesfull")}
}

plotNETWORK <- function(gg, mem, out_dir, fname, match, p, co_ind) {

  png(filename = paste0(out_dir, "/", fname, "_peak", p, "_and_CoPeaks.png"),
      width = 10, height = 8, units = "in", res = 100)

  ## PID of the peak of interest
  pmain <- paste0("Main peak: ", round(pks[p, "mz"], digits = 4),
                  ". Coeluting peaks: ",
                  length(co_ind))
  psub <- paste0("Features in POI component: ",
                 length(which(mem == 1 | mem == 2)))

  colors <- c("#ABDDA4", "#E6F598", "grey")

  ## this could be used for alternating colors of edges
  # edge_colors <- c("royalblue4" , "grey")
  # edg <- ifelse(igraph::E(gg)$weight > thr, 1, 2)

  ## make coordinates for vertices based on correlation: scale cor coef to maximise distance
  ## using Fruchterman-Reingold layout algorithm to prevent overlap
  coords <- igraph::layout_(gg, igraph::with_fr(weights = scaleEDGES(igraph::E(gg)$weight, from = 0.01, to = 10)))

  igraph::plot.igraph(gg,
                      main = pmain,
                      sub = psub,
                      layout = coords,
                      vertex.color = colors[mem],
                      vertex.size = 22,
                      vertex.label.color = "black",
                      # edge.color = edge_colors[edg],
                      edge.color = "royalblue4",
                      edge.label = round(igraph::E(gg)$weight, digits = 2),
                      # edge.label.color = edge_colors[edg],
                      edge.label.color = "royalblue4")
  dev.off()
}

buildCOR <- function(co_ind, eic, pearson) {
  poi_co_cor <- expand.grid(x = co_ind, y = co_ind) %>%
    group_by(x, y) %>%
    mutate(cor = massflowR::getCOR(x = x, y = y, eic = eic, pearson = pearson)) %>%
    filter(x != y)
}

scaleEDGES <- function(x, from, to) {
  (x - min(x)) / max(x - min(x)) * (to - from) + from
  }


