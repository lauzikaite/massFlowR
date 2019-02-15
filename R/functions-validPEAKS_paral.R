####---- function for bplapply call, not to export
validPEAKS_paral <-
  function(pkg,
           peakgrs_ints,
           out_dir,
           cor_thr = 0.75,
           min_samples_n = 3,
           save_plot = FALSE
           ) {
    ## extract intensities for the peak-group of interest
    pkg_ints <- peakgrs_ints[[pkg]]
    
    ####---- (A) build correlation network between all peaks in the peak-group ----
    peaks_ids <- unique(pkg_ints$peakid)
    ## make peak-peak pairs
    pkg_cor <- setNames(as.data.frame(t(utils::combn(
      peaks_ids, 2, simplify = T
    ))),
    ## column names "from" and "to" will be needed in graph generation
    nm = c("from", "to"))
  
    ## weight is the correlation coefficient
    pkg_cor$weight <- apply(
      pkg_cor,
      1,
      FUN = corPEAKS,
      pkg_ints = pkg_ints,
      min_samples_n = min_samples_n
    )
    
    ## build network between peaks in the peak-group and report communities
    title <- paste0("Peak-group-", pkg)
    pkg_coms <-
      buildGRAPH(
        pkg_cor = pkg_cor,
        cor_thr = cor_thr,
        title = title,
        out_dir = out_dir,
        plot = save_plot
      )
    
    ## split peaks into new peak-groups according to communities
    pkg_ints$new_peakgr <-
      sapply(pkg_ints$peakid, function(peakid) {
        pkg_coms[[match(peakid, names(pkg_coms))]]
      })

    # ## plot intensities across all samples
    # colors <-
    #   setNames(viridis::magma(
    #     end = 0.9,
    #     ## avoid bright yellow
    #     begin = 0.2,
    #     ## avoid dark purple
    #     n = length(unique(pkg_ints$new_peakgr))
    #   ),
    #   nm = unique(pkg_ints$new_peakgr))
    # 
    # ## linetype according to how many peaks in a group
    # peaks_per <- sapply(pkg_ints$new_peakgr, function(peakgr) {
    #   length(unique(pkg_ints[which(pkg_ints$new_peakgr == peakgr), "peakid"]))
    # })
    # linetype <- setNames(sapply(peaks_per, function(p) {
    #   if (p > 2) {
    #     "solid"
    #   } else {
    #     if (p == 2) {
    #       "dashed"
    #     } else {
    #       "dotted"
    #     }
    #   }
    # }), nm = sapply(peaks_per, function(p) {
    #   if (p > 2) {
    #     ">2"
    #   } else {
    #     if (p == 2) {
    #       2
    #     } else {
    #       1
    #     }
    #   }
    # }))
    # pkg_ints$peaks_per <- names(linetype)
    # 
    # if (save_plot == TRUE) {
    #   gg <-  ggplot2::ggplot(data = pkg_ints) +
    #     ggplot2::facet_wrap(~ new_peakgr) +
    #     ggplot2::geom_line(
    #       ggplot2::aes(
    #         x = run_order,
    #         y = into,
    #         linetype = as.factor(peaks_per),
    #         color = as.factor(new_peakgr),
    #         group = as.factor(peakid)
    #       )
    #     ) +
    #     ggplot2::scale_color_manual(values = colors,
    #                                 name = "Peaks sub-groups ") +
    #     ggplot2::scale_linetype_manual(name = "Unique peaks per group",
    #                                    values = linetype) +
    #     ggplot2::ylab("Intensity") +
    #     ggplot2::xlab("Sample run order") +
    #     ggplot2::ggtitle(title) +
    #     ggplot2::theme_bw()
    # 
    #   grid::grid.newpage()
    #   grDevices::pdf(
    #     width = 8,
    #     height = 8,
    #     file = paste0(
    #       out_dir,
    #       "/",
    #       "peak-group-",
    #       pkg,
    #       "_intensity-across-samples.pdf"
    #     )
    #   )
    #   grid::grid.draw(ggplot2::ggplotGrob(gg))
    #   grDevices::dev.off()
    # }
    
    ####---- (B) extract intensity data for each detected community and return as a list of data frames ----
    # find best-representative peaks in all split-peak-groups
    # peak_counts <- lapply(unique(pkg_ints$new_peakgr), FUN = countPEAKS, pkg_ints = pkg_ints)
    new_peakgrs <- lapply(unique(pkg_ints$new_peakgr), function(npg) {
      new_pkg <- pkg_ints[pkg_ints$new_peakgr == npg,]
      # drop temporaly assigned new_peakgr id, since these are not unique across all peak-groups
      new_pkg[, !(names(new_pkg) %in% "new_peakgr")]
    })
    
    return(new_peakgrs)
    
  }



####---- extract peak-group data from every sample
extractPEAKGR <- function(pkg, object, samples) {
  pkg_peakids <- object@tmp$peakid[which(object@tmp$peakgr == pkg)]

  ## for every sample, extract mz/rt/into info for  each peak in the peak-group
  snames <- names(object@data)
  pkg_ints <- lapply(snames, function(s) {
    sample <- object@data[[s]]
    sample <-
      setNames(sample[which(sample$tmp_peakid %in% pkg_peakids), c("tmp_peakid", "into", "mz", "rt")],
               c("peakid", "into", "mz", "rt"))
    sample$run_order <-
      rep(object@samples[which(snames == s), "run_order"], nrow(sample))

    return(sample)
  })
  pkg_ints <- do.call("rbindCLEAN", pkg_ints)
  pkg_ints$peakgr <- rep(pkg, nrow(pkg_ints))
  return(pkg_ints)

}



####---- extract peak-group intensities from every sample
# extractINT <- function(pkg, object, samples) {
#   pkg_peakids <- object@tmp$peakid[which(object@tmp$peakgr == pkg)]
#   ## for every sample, extract intensities for peakgr of interest
#   snames <- names(object@data)
#   pkg_ints <- lapply(snames, function(s) {
#     sample <- object@data[[s]]
#     sample <-
#       setNames(sample[which(sample$tmp_peakid %in% pkg_peakids), c("tmp_peakid", "into")],
#                c("peakid", "into"))
#     sample$filepath <- rep(s, nrow(sample))
#     sample$run_order <-
#       rep(object@samples[which(snames == s), "run_order"],
#           nrow(sample))
#     sample$sample_no <- rep(which(snames == s), nrow(sample))
#     return(sample)
#   })
# 
#   pkg_ints <- do.call("rbindCLEAN", pkg_ints)
#   pkg_ints$peakgr <- rep(pkg, nrow(pkg_ints))
#   return(pkg_ints)
# }

####---- extract peak-group intensities from every sample
## add NA for peaks that are missing in a sample
# extractINTv2 <- function(pkg, object) {
#   peakid <- object@tmp$peakid[which(object@tmp$peakgr == pkg)]
#   snames <- names(object@data)
#   
#   pkg_ints <- lapply(snames, function(s) {
#     sample <- object@data[[s]]
#     sample_order <- object@samples[which(snames == s), "run_order"]
#     sample_pkg <- sample[which(sample$tmp_peakgr == pkg), ]
#     sample_out <- data.frame(
#       peakid = peakid,
#       filepath = rep(s, length(peakid)),
#       run_order = rep(sample_order, length(peakid)),
#       into = NA,
#       stringsAsFactors = F
#     )
#     
#     ## add into values from sample
#     common_tmp <-
#       match(sample_pkg$tmp_peakid, sample_out$peakid, nomatch = FALSE)
#     sample_out[common_tmp, "into"] <- sample_pkg$into
#     return(sample_out)
#   })
#   
#   pkg_ints <- do.call("rbindCLEAN", pkg_ints)
#   pkg_ints$peakgr <- rep(pkg, nrow(pkg_ints))
#   return(pkg_ints)
# }






####---- correlate two peaks' intensities across all samples
corPEAKS <- function(pair, pkg_ints, min_samples_n) {
  x <- pair["from"]
  y <- pair["to"]
  
  pkg_ints_x <- pkg_ints[which(pkg_ints$peakid == x),]
  pkg_ints_y <- pkg_ints[which(pkg_ints$peakid == y),]
  common_samples <-
    base::intersect(pkg_ints_x$run_order, pkg_ints_y$run_order)
  
  ## use min_samples parameter to omit peaks that are not present in enough samples
  if (length(common_samples) > min_samples_n) {
    ## if multiple peaks in a single sample were matched to the same tmp peak
    ## correlate only the most intense
    ix <- sapply(common_samples,
                 FUN = getINT,
                 peak = pkg_ints_x)
    
    iy <- sapply(common_samples,
                 FUN = getINT,
                 peak = pkg_ints_y)
    cc <- cor(ix, iy, method = "pearson")
    ## negative coeficients would break graph generation
    cc <- ifelse(cc < 0, 0, cc)
  } else {
    cc <- 0
  }
  return(cc)
}

getINT <- function(s, peak) {
  ints <- peak$into[which(peak$run_order == s)]
  ints <- ints[which(ints == max(ints))]
  return(ints)
}



# getPCS <- function(n, peakgrs_split_sel) {
#   
#   pkg <- peakgrs_split_sel[[n]]
#   
#   ## get mz&rt averages for eack peak in every peak-group split
#   peaks <- unique(pkg$peakid)
#   pkg_out <- lapply(peaks, function(p) {
#     peak <- pkg[pkg$peakid == p, ]
#     peak_out <- data.frame(peakid = p,
#                            mz = mean(peak$mz),
#                            rt = mean(peak$rt),
#                            peakgr = unique(peak$peakgr))
#     return(peak_out)
#   })
#   pkg_out <- do.call("rbindCLEAN", pkg_out)
#   pkg_out$pcs <- n
#   return(pkg_out)
# }


extractPCS <- function(sname, object, final_tmp) {
  
  ## take original peak picking data of the sample
  snames <- object@samples$filename
  sn <- which(snames == sname)
  sdata <- object@data[[sn]]

  ## get all of the peaks from the final list
  peaks <- match(final_tmp$peakid, sdata$tmp_peakid)
  sdata_out <- as.data.frame(sdata[peaks,"into"],
                             row.names = NULL)
  colnames(sdata_out) <- as.character(sname)
  return(sdata_out)
}

# extractPeaks <- function(peakid, object) {
#   
#   ## get mz, rt, mzmin, mzmax, rtmin, rtmax values from every sample
#   for (sn in 1:nrow(samples)) {
#     
#     sdata <- object@data[[sn]]
#     sdata[match(peakid, sdata$tmp_peakid), c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax"  )]
#     
#   }
#   
#   
# }
