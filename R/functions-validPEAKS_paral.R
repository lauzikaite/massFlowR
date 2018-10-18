####---- function for bplapply call, not to export
validPEAKS_paral <- function(pkg, peakgrs_ints, out_dir, cor_thr = 0.75, min_samples = 10) {
  
  ## extract intensities for the peak-group of interest
  pkg_ints <- peakgrs_ints[[pkg]]
  
  ## build correlation matrix between all peaks in the peak-group
  peaks_ids <- unique(pkg_ints$peakid) ## extract peak-ids for the peak-group
  pkg_cor <- setNames(as.data.frame(t(combn(peaks_ids, 2, simplify = T))), nm = c("from", "to"))  ## make peak-peak pairs
  min_samples_n <- ceiling(length(unique(pkg_ints$run_order)) * (min_samples/100))
  pkg_cor$weight <- apply(pkg_cor, 1, FUN = corPEAKS, pkg_ints = pkg_ints, min_samples_n = min_samples_n) ## weight is the correlation coefficient

  ## build network between peaks in the peak-group and detect communities
  pkg_coms <- buildGRAPH(cor_dt = pkg_cor, cor_thr = cor_thr)
  
  ## split peaks into new peak-groups according to communities
  pkg_ints$new_peakgr <- sapply(pkg_ints$peakid, function(peakid) {
    pkg_coms[[match(peakid, names(pkg_coms))]]
  })
    
  ## plot intensities across all samples
  colors <- setNames(viridis::viridis(end = 0.9, ## avoid bright yellow
                                      begin = 0.2, ## avoid dark purple
                                      n = length(unique(pkg_ints$new_peakgr))),
                     nm = unique(pkg_ints$new_peakgr))
  gg <-  ggplot2::ggplot(data = pkg_ints) +
    ggplot2::geom_line( ggplot2::aes(x = run_order,
                                     y = into,
                                     color = as.factor(new_peakgr),
                                     group = as.factor(peakid))) +
    ggplot2::scale_color_manual(values = colors,
                                name = "Peaks sub-groups ") +
    ggplot2::ylab("Intensity") +
    ggplot2::xlab("Sample run order") +
    ggplot2::theme_bw()
  
  grid::grid.newpage()
  grDevices::pdf(width = 8, height = 8,
                 file = paste0(out_dir, "/", "peak-group-", pkg, "_intensity-across-samples.pdf"))
  grid::grid.draw(ggplot2::ggplotGrob(gg))
  grDevices::dev.off()
  
  return(pkg_ints)

}
 
  


## extract peak-group intensities from every sample
extractINT <- function(pkg, object) {
  pkg_peakids <- object@tmp$peakid[which(object@tmp$peakgr == pkg)]
  pkg_ints <- lapply(names(object@data), function(s) {
    sample <- object@data[[s]]
    sample <- sample[which(sample$peakid %in% pkg_peakids),c("peakid","into")]
    sample$filepath <- rep(s, nrow(sample))
    sample$run_order <- rep(object@samples[which(object@samples == s),"run_order"], nrow(sample))
    return(sample)
  })
  pkg_ints <- do.call(function(...) rbind(..., make.row.names = F), pkg_ints)
  pkg_ints$peakgr <- rep(pkg, nrow(pkg_ints))
  return(pkg_ints)
}

## Correlate two peaks' intensities across all samples
corPEAKS <- function(pair, pkg_ints, min_samples_n) {
  x <- pair["from"]
  y <- pair["to"]
  sx <- pkg_ints[which(pkg_ints$peakid == x),]
  sy <- pkg_ints[which(pkg_ints$peakid == y),]
  common_samples <- base::intersect(sx$filepath, sy$filepath)
  if (length(common_samples) > min_samples_n) {
    ix <- sx$into[which(sx$filepath %in% common_samples)]
    iy <- sy$into[which(sy$filepath %in% common_samples)]
    cc <- cor(ix, iy, method = "pearson")
    cc <- ifelse(cc < 0, 0, cc) ## negative coeficients would break graph generation
  } else {
    cc <- 0
  }
  return(cc)
}

## build graph
buildGRAPH <- function(cor_dt, cor_thr) {
  
  g <- igraph::graph_from_data_frame(cor_dt, directed = FALSE)
  
  ## delete edges between peak pairs with cor < thr
  g <- igraph::delete.edges(g, igraph::E(g)[[which(igraph::E(g)$weight < cor_thr)]])
  
  ## detect community structure in the network where non-correlated pairs are omitted
  coms <- igraph::cluster_label_prop(g)
  mem <- igraph::membership(coms)
  return(mem)
}
