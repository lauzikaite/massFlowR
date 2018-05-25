#' Title
#'
#' @param Pearson
#' @param match
#' @param thr
#' @param plot
#' @param pks
#'
#' @return
#' @export
#'
#' @examples
buildCOMPS <- function(Pearson, match, thr = 0.95, plot = FALSE, pks, clean = TRUE) {

  print(paste0("Apex match by scpos window: ", match))
  print(paste0("Pearson correlation: ", Pearson))

  start <- Sys.time()

  ## make a duplicate reduced peaks table which will be used in building of the final output df all_poi
  pks_sel <- pks %>%
    select(pid, mz:maxo, poi_name)

  ## make a duplicate reduced peaks table, where componenents will be stored for every peak
  pkscomps <- pks %>%
    mutate(comp = NA) %>%
    select(pid = pid, mz, scpos, into, comp) %>%
    data.frame

  # all_poi <- data.frame()

  ####--- (A) while there is still a peak that was not assigned to component yet ----
  # while(any(is.na(pkscomps$comp))) { p <- pkscomps %>% filter(is.na(comp)) %>% slice(1) %>% pull(pid)

  pids <- pkscomps %>%
    pull(pid)

  for (p in pids) {
    if (is.na(pkscomps[p, "comp"])) {

      ####---- (A1) find all co-eluting peaks which are not asigned to a component yet and match by their scpos being in the POI scpos ----
      if (match == 0) {
        ## version(1) - exact scpos match
        pkscomps_sel <- pkscomps %>%
          filter(is.na(comp) & scpos == scpos[pid == p])

      } else {
        ## version(2) - scpos +- match
        mat <- c(pkscomps[p, "scpos"] - match, pkscomps[p, "scpos"] + match)
        pkscomps_sel <- pkscomps %>%
          filter(is.na(comp)) %>% filter(between(scpos, mat[1], mat[2]))
      }

      co_ind <- pkscomps_sel %>%
        pull(pid)

      ####--- if POI has any co-eluting peaks in the peak table ----
      if (length(co_ind) > 1 & p %in% co_ind) {

        print(paste0("Checking p: ", p))

        ####--- build correlation matrix for all co-eluting peaks ----

        ## build a cor mat between all co-eluting feature pairs, correlate each pair separately by leaving only common scans
        poi_co_cor <- expand.grid(x = co_ind, y = co_ind) %>%
          group_by(x, y) %>%
          mutate(cor = getCOR(x = x, y = y)) %>%
          filter(x != y)

        ####---- build interaction network with full co-eluting set ----

        ## if there's at least one CO pair with cor above threshold
        if (any(poi_co_cor$cor > thr)) {

          pkscomps <- buildNETWORK(poi_co_cor = poi_co_cor, pkscomps = pkscomps, co_ind = co_ind, thr = thr, pks = pks, p = p, plot = plot)

          } else next
        } else next
      } else next
  }

  ## assign component id to every remaining peak, which was not correlated to anything else
  pkscomps <- merge(pks, pkscomps[, c("pid", "comp")],
                    by = c("pid"), all = T)
  pkscomps[which(is.na(pkscomps$comp)), "comp"] <- seq(
    from = (max(na.omit(pkscomps$comp)) + 1),
    length.out = length(which(is.na(pkscomps$comp))),
    by = 1)

  if (clean == TRUE) {
    message("'clean' set to TRUE. Removing 1-peak components")
    pkscomps <- removeCOMPS(pkscomps = pkscomps)
  }

  print(Sys.time() - start)
  return(pkscomps)
}


####---- helper functions, not to export ----
removeCOMPS <- function(pkscomps) {

  comps1 <- pkscomps %>%
    group_by(comp) %>%
    summarise(n = n()) %>%
    filter(n == 1) %>%
    pull(comp)

  pkscomps <- pkscomps %>%
    filter(!comp %in% comps1) %>%
    mutate(order_ID = row_number())

  return(pkscomps)

}

getCOR <- function(x, y) {

  rx <- MSnbase::rtime(eic[[x]])
  ry <- MSnbase::rtime(eic[[y]])
  common_scan <- base::intersect(rx, ry)
  if (length(common_scan) > 3) {
    ix <- as.numeric(MSnbase::intensity(eic[[x]])[which(rx %in% common_scan)])
    iy <- as.numeric(MSnbase::intensity(eic[[y]])[which(ry %in% common_scan)])
    cc <- ifelse(Pearson == FALSE, cor(ix, iy, method = "spearman"), cor(ix, iy, method = "pearson"))
  } else {
    cc <- 0
  }
  return(cc)
}

plotNETWORK <- function(gg, mem, out_dir_fname, fname, p) {

  png(filename =
        paste0(out_dir_fname, "/", fname, "_pks", p, "_and_CoPeaks.png"),
      width = 10, height = 8, units = "in", res = 100)

  ## PID of the peak of interest
  pmain <- paste0("Main peak: ", round(pks[p, "mz"], digits = 4),
                  ". Coeluting peaks: ",
                  length(co_ind))
  psub <- paste0("Features in POI component: ",
                 length(which(mem == 1 | mem == 2)))

  colors <- c("#9cc057", "#cad587", "grey")
  igraph::plot.igraph(gg,
                      main = pmain,
                      sub = psub,
                      layout = layout.fruchterman.reingold,
                      vertex.color = colors[mem],
                      vertex.size = 22,
                      vertex.label.color = "black",
                      edge.color = "black",
                      edge.width = igraph::E(gg)$weight)
  dev.off()
}

buildNETWORK <- function(poi_co_cor, pkscomps, co_ind, thr, pks, p, plot, out_dir_fname, fname) {

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
    plotNETWORK(gg = gg, mem = mem, out_dir_fname = out_dir_fname, fname = fname, p = p)
  }

  ####---- if there is at least once community with more than one member, asign community's compounds to a component ----
  members <- which(mem == 2 | mem == 1)

  if (length(members) >= 2) {

    component <- co_ind[members]

    ####---- save components metadata: feature pairs, their cor and peak info ----
    # gg <- igraph::delete.vertices(gg, c(which(mem == 3)))
    #
    # pairs <- setNames(
    #   data.frame(igraph::as_edgelist(gg, names = F)), c("x", "y")) %>%
    #   mutate(
    #     xind = co_ind[members[x]], yind = co_ind[members[y]],
    #     pair_ID = paste(xind, yind, sep = "_"))
    # co_pairs <- poi_co_cor %>%
    #   filter(cor >= thr) %>%
    #   filter(paste(x, y, sep = "_") %in% pairs$pair_ID) %>% ## add component number to each feature: redundant, but helps to spot if component asignment is correct
    #   left_join(
    #     pkscomps %>%
    #       select(pid, comp_x = comp),
    #     by = c(x = "pid")) %>%
    #   left_join(
    #     pkscomps %>%
    #       select(pid, comp_y = comp),
    #     by = c(y = "pid")) %>%
    #   ## add remaining peaks information, e.g. mz, rt ...
    #   left_join(pks_sel, by = c(x = "pid")) %>%
    #   left_join(pks_sel, by = c(y = "pid"), suffix = c("_x", "_y"))

    # all_poi <- co_pairs %>% data.frame %>% rbind(all_poi)

    ####---- update table with assigned component IDs ----
    comp <- ifelse(all(is.na(pkscomps$comp)), 1, max(pkscomps$comp, na.rm = T) + 1)
    pkscomps[component, "comp"] <- rep(comp, length(component))
  }
  return(pkscomps)
}
