####---- function for bplapply call, not to export
validPEAKS_paral <- function(pkg, peakgrs_ints, out_dir, cor_thr = 0.75, min_samples = 10) {
  
  ## extract intensities for the peak-group of interest
  pkg_ints <- peakgrs_ints[[pkg]]
  
  ## build correlation matrix between all peaks in the peak-group
  peaks_ids <- unique(pkg_ints$peakid) ## extract peak-ids for the peak-group
  pkg_cor <- setNames(as.data.frame(t(combn(peaks_ids, 2, simplify = T))),  ## make peak-peak pairs
                      nm = c("from", "to")) ## column names "from" and "to" will be needed in graph generation
  
  min_samples_n <- ceiling(length(unique(pkg_ints$run_order)) * (min_samples/100))
  pkg_cor$weight <- apply(pkg_cor, 1,
                          FUN = corPEAKS,
                          pkg_ints = pkg_ints, min_samples_n = min_samples_n) ## weight is the correlation coefficient

  ## build network between peaks in the peak-group and detect communities
  title <- paste("Peak-group no:", pkg)
  pkg_coms <- buildGRAPH(pkg_cor = pkg_cor, cor_thr = cor_thr, title = title, out_dir = out_dir, pkg = pkg)
  
  ## split peaks into new peak-groups according to communities
  pkg_ints$new_peakgr <- sapply(pkg_ints$peakid, function(peakid) {
    pkg_coms[[match(peakid, names(pkg_coms))]]
  })
    
  ## plot intensities across all samples
  colors <- setNames(viridis::magma(end = 0.9, ## avoid bright yellow
                                      begin = 0.2, ## avoid dark purple
                                      n = length(unique(pkg_ints$new_peakgr))),
                     nm = unique(pkg_ints$new_peakgr))
  
  ## linetype according to how many peaks in a group
  peaks_per <- sapply(pkg_ints$new_peakgr, function(peakgr) {
    length(unique(pkg_ints[which(pkg_ints$new_peakgr == peakgr), "peakid"]))
  })
  linetype <- setNames(sapply(peaks_per, function(p) {
    if (p > 2) {
      "solid"
    } else {
      if (p == 2) {
        "dashed"
      } else {
        "dotted"
      }
    }
  }), nm = sapply(peaks_per, function(p) {
    if (p > 2) {
      ">2"
    } else {
      if (p == 2) {
        2
      } else {
        1
      }
    }
  }))
  pkg_ints$peaks_per <- names(linetype)
  gg <-  ggplot2::ggplot(data = pkg_ints) +
    ggplot2::geom_line( ggplot2::aes(x = run_order,
                                     y = into,
                                     linetype = as.factor(peaks_per),
                                     color = as.factor(new_peakgr),
                                     group = as.factor(peakid))) +
    ggplot2::scale_color_manual(values = colors,
                                name = "Peaks sub-groups ") +
    ggplot2::scale_linetype_manual(name = "Unique peaks per group",
                                   values = linetype) +
    ggplot2::ylab("Intensity") +
    ggplot2::xlab("Sample run order") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw()

  grid::grid.newpage()
  grDevices::pdf(width = 8, height = 8,
                 file = paste0(out_dir, "/", "peak-group-", pkg, "_intensity-across-samples.pdf"))
  grid::grid.draw(ggplot2::ggplotGrob(gg))
  grDevices::dev.off()
  
  return(pkg_ints)

}
 
  


####---- extract peak-group intensities from every sample
extractINT <- function(pkg, object) {
  pkg_peakids <- object@tmp$peakid[which(object@tmp$peakgr == pkg)]
  pkg_ints <- lapply(names(object@data), function(s) {
    sample <- object@data[[s]]
    sample <- setNames(sample[which(sample$tmp_peakid %in% pkg_peakids),c("tmp_peakid","into")], c("peakid", "into"))
    sample$filepath <- rep(s, nrow(sample))
    sample$run_order <- rep(object@samples[which(object@samples == s),"run_order"], nrow(sample))
    return(sample)
  })
  pkg_ints <- do.call(function(...) rbind(..., make.row.names = F), pkg_ints)
  pkg_ints$peakgr <- rep(pkg, nrow(pkg_ints))
  return(pkg_ints)
}

####---- correlate two peaks' intensities across all samples
corPEAKS <- function(pair, pkg_ints, min_samples_n) {
  x <- pair["from"]
  y <- pair["to"]
  
  x <- pkg_ints[which(pkg_ints$peakid == x),]
  y <- pkg_ints[which(pkg_ints$peakid == y),]
  common_samples <- base::intersect(x$filepath, y$filepath)
  
  ## use min_samples parameter to omit peaks that are not present in enough samples
  if (length(common_samples) > min_samples_n) {
    
    ## if multiple peaks in a single sample were matched to the same tmp peak
    ## correlate only the most intense
    ix <- sapply(common_samples,
                 FUN = getINT,
                 peak = x)
    
    iy <- sapply(common_samples,
                 FUN = getINT,
                 peak = y)
    cc <- cor(ix, iy, method = "pearson")
    cc <- ifelse(cc < 0, 0, cc) ## negative coeficients would break graph generation
  } else {
    cc <- 0
  }
  return(cc)
}


getINT <- function(s, peak) {
  ints <- peak$into[which(peak$filepath == s)]
  ints <- ints[which(ints == max(ints))]
  return(ints)
}



#####---- build graph of peaks' intensities correlation across all samples
buildGRAPH <- function(pkg_cor, cor_thr, title, out_dir, pkg) {
  
  g <- igraph::graph_from_data_frame(pkg_cor, directed = FALSE)

  # colors <- c("royalblue4", "lightgrey", "grey")
  # edge_colors <- sapply(igraph::E(g)$weight, function(w) {
  #   
  #   if (w > cor_thr) {
  #     colors[1]
  #   } else {
  #     if (w > 0.5) {
  #       colors[2]
  #     } else {
  #     colors[3]
  #     }
  #   }
  # })
  # igraph::E(g)[1:(length(E(g)))]$color <-  edge_colors
  
  ## scale edge thickness 
  igraph::E(g)$weight <- scaleEDGES(igraph::E(g)$weight, from = 0.01, to = 10)
  
  ## make coordinates for vertices based on correlation: scale cor coef to maximise distance
  ## using Fruchterman-Reingold layout algorithm to prevent overlap
  coords <- igraph::layout_with_fr(g)
  
  ## delete edges between peak pairs with cor < thr
  g <- igraph::delete.edges(g, igraph::E(g)[[which(igraph::E(g)$weight < cor_thr)]])
  
  ## detect community structure in the network where non-correlated pairs are omitted
  coms <- igraph::cluster_label_prop(g)
  mem <- igraph::membership(coms)
  
  grid::grid.newpage()
  grDevices::pdf(width = 8, height = 8,
                 file = paste0(out_dir, "/", "peak-group-", pkg, "_network.pdf"))
  igraph::plot.igraph(g,
                      main = title,
                      sub = paste("Cor threshold:", cor_thr),
                      layout = coords)
  grDevices::dev.off()
  
  return(mem)
}
