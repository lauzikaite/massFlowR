# extractPEAKGR ---------------------------------------------------------------------------------------------------
#' @title Extract peak-group data from every sample
#'
#' @param pkg \code{numeric} specifying a peak-group-of-interest.
#' @param object \code{massFlowTemplate} class object.
#' @param samples \code{data.frame} with experiment's metadata.
#'
#' @return Function returns a \code{data.frame} with centWave values for the peak-group-of-interest in every sample. 
#' 
extractPEAKGR <- function(pkg, object, samples) {
  pkg_peakids <- object@tmp$peakid[which(object@tmp$peakgr == pkg)]
  
  ## for every sample, extract mz/rt/into info for  each peak in the peak-group
  snames <- names(object@data)
  req_cnames <- c("tmp_peakid", "into", "maxo", "mz", "rt")
  new_cnames <- c("peakid", "into", "maxo", "mz", "rt")
  pkg_ints <- lapply(snames, function(s) {
    sample <- object@data[[s]]
    sample <-
      setNames(sample[which(sample$tmp_peakid %in% pkg_peakids), req_cnames], nm = new_cnames)
    sample$run_order <-
      rep(object@samples[which(snames == s), "run_order"], nrow(sample))
    return(sample)
  })
  pkg_ints <- do.call("rbindCLEAN", pkg_ints)
  pkg_ints$peakgr <- rep(pkg, nrow(pkg_ints))
  return(pkg_ints)
}


# validPEAKGR ---------------------------------------------------------------------------------------------------
#' @title Validate a peak-group-of-interest by performing a cross-sample Pearson correlation.
#' 
#' @description Functions correlates all peaks of a peak-groups across all samples in which they were detected.
#' Obtained Pearson correlation coefficients are used to build a peak-peak network with pairs above selected correlation threshold.
#' \emph{igraph} Cluster Label Propagation algorithm is used to detect communities of peaks with the network.
#' These communities are returned as independent peak-goups, or, Pseudo Chemical Spectra.
#'
#' @param pkg \code{numeric} specifying a peak-group-of-interest.
#' @param pkg_ints  \code{data.frame} with centWave values for the peak-group-of-interest.
#' @param out_dir \code{character} specifying desired directory for output, if plotting is selected.
#' @param cor_thr \code{numeric} defining Pearson correlation coefficient threshold for inter-sample correlation between peaks. 
#' @param min_samples_n \code{numeric} specifying the minimum number of samples in which peak has to be detected in order to be considered.
#' @param save_plot \code{logical} whether network plot should be saved as a pdf file in the out_dir, default set to FALSE.
#'
#' @return Function returns a \code{list} for every community into which peak-group was split to by the Cluster Label Propagation algorithm.
#' Each element of the returned list contains centWave values for every community peak in the samples in which it was detected.
#'
validPEAKGR <-
  function(pkg,
           pkg_ints,
           out_dir,
           cor_thr,
           min_samples_n,
           save_plot = FALSE
           ) {
    ####---- build correlation network between all peaks in the peak-group ----
    pkg_cor <- t(utils::combn(unique(pkg_ints$peakid),
                              2,
                              simplify = T))
    
    ## weight is the correlation coefficient
    weight <- apply(
      pkg_cor,
      1,
      FUN = corPEAKS,
      pkg_ints = pkg_ints,
      min_samples_n = min_samples_n
    )
    
    ## build network between peaks in the peak-group and report communities
    ## data.frame is required for igraph in buildGRAPH
    pkg_cor <-
      setNames(data.frame(cbind(pkg_cor, weight)), nm = c("from", "to", "weight"))
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

    ####---- extract intensity data for each detected community and return as a list of data frames ----
    new_peakgrs <- lapply(unique(pkg_ints$new_peakgr), function(npg) {
      new_pkg <- pkg_ints[pkg_ints$new_peakgr == npg,]
      # drop temporaly assigned new_peakgr id, since these are not unique across all peak-groups
      new_pkg[, !(names(new_pkg) %in% "new_peakgr")]
    })
    return(new_peakgrs)
  }


# corPEAKS --------------------------------------------------------------------------------------------------------
#' @title Correlate two peaks' intensities across all samples.
#' 
#' @description Function takes a pair of peaks and performs Pearson correlation of their "into" values across the samples in which they were detected.
#' 
#' @param pair \code{matrix} with 1 row and 2 columns: "from" and "to", which indicate peaks' peakids.
#' @param pkg_ints \code{data.frame} which must contain columns "peakid", "run_order" and "into", 
#' indicating intensities of every peak in a peak-group in every sample in which they were detected.
#' @param min_samples_n \code{numeric} specifying the minimum number of samples in which peak has to be detected in order to be correlated with other peaks.
#'
#' @return Function returns a \code{numeric} indicating Pearson correlation coefficient between the two peaks of interest.
#' 
corPEAKS <- function(pair, pkg_ints, min_samples_n) {
  x <- pair[1]
  y <- pair[2]
  
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
    cc <- cor(ix, iy, method = "pearson", use = "pairwise.complete.obs")
    ## negative coeficients would break graph generation
    cc <- ifelse(cc < 0, 0, cc)
  } else {
    cc <- 0
  }
  return(cc)
}


# getINT --------------------------------------------------------------------------------------------------------
#' @title Get intensity of a peak-of-interest in a single sample.
#' 
#' @description Get intensity of a peak-of-interest in a single sample.
#'
#' @param s \code{numeric} indicating sample's run order.
#' @param peak \code{data.frame} with intensities of the peak-of-interest in all samples.
#'
#' @return Function returns an \code{integer} with peak intensities.
#' 
getINT <- function(s, peak) {
  ints <- peak$into[which(peak$run_order == s)]
  ints <- ints[which(ints == max(ints))]
  return(ints)
}


# extractCOMMUNITIES ----------------------------------------------------------------------------------------------
#' @title Extract communities, or, Pseudo Chemical Spectra, with > 1 peak
#' 
#' @description Function extracts communities, or, Pseudo Chemical Spectra, with > 1 peak, and returns final peak table.
#'
#' @param peakgrs_split \code{list} for every community into which peak-group was split to by the Cluster Label Propagation algorithm.
#'
#' @return Function returns a \code{data.frame} with peakids and corresponding Pseudo Chemical Spectra (PCS).
#' 
#' @seealso \code{\link{validPEAKGR}}, \code{\link{validPEAKS}}
#' 
extractCOMMUNITIES <- function(peakgrs_split) {
  ## iterate over every peak-group
  peakgrs_comms_peakids <- lapply(peakgrs_split, FUN = getPEAKIDS)
  
  ## retain only communities with >1 peaks
  comms <-
    Filter(length, unlist(peakgrs_comms_peakids, recursive = F))
  
  ## convert list to table
  comms_dt <- data.frame(peakid = unlist(comms),
                         pcs = rep(1:length(comms), sapply(comms, length)))
  return(comms_dt)
}


# getPEAKIDS ------------------------------------------------------------------------------------------------------
#' @title Extract peakids from every community/Pseudo Chemical Spectra of a peak-group-of-interest.
#' 
#' @description Function is used after \code{\link{validPEAKGR}} has been applied to identify communities of peaks using correlation network analysis.
#' Function returns peakids of every community into which a peak-group-of-interest was split to.
#'
#' @param pkg \code{list} with communities into which each peak-group was split to.
#' List length is equal to number of peak-groups.
#'
#' @return Function returns peakids for the given network community.
#' 
#' @seealso \code{\link{validPEAKGR}}, \code{\link{validPEAKS}}.
#' 
getPEAKIDS <- function(pkg) {
  ## iterate over every community 
  pkg_comms <- lapply(pkg, function(comm) {
    peaks <- unique(comm$peakid)
    if (length(peaks) > 1) {
      return(peaks)
    } else {
      NULL
    }
  }
  )
  pkg_comms_sel <- pkg_comms[which(sapply(pkg_comms, length) > 0)]
  return(pkg_comms_sel)
}



# exportSAMPLE ----------------------------------------------------------------------------------------------
#' @title Extract centWave values for peaks in the final template for a sample-of-interest.
#' 
#' @description Function extracts centWave peak-picking values for a single sample.
#' Only peaks in the final template are returned.
#' 
#' @param sdata \code{data.frame} with centWave measures for the sample-of-interest.
#' @param final_tmp \code{data.frame} with columns "peakid" and "pcs" indicating to which Pseudo Chemical Spectra each peak (peakid) belongs to.
#' Contains only the peaks that were returned by peak validation algorithm.
#'
#' @return Function returns a \code{data.frame} with all centWave values for all peaks from the final template that were detected in the sample-of-interest.
#' 
#' @seealso \code{\link{validPEAKGR}}, \code{\link{validPEAKS}}
#' 
exportSAMPLE <- function(sdata, final_tmp) {
  ## get all of the peaks from the final list
  peaks <- match(final_tmp$peakid, sdata$tmp_peakid)
  sdata_out <- cbind(peakid = final_tmp$peakid,
                     sdata[peaks, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo")])
  rownames(sdata_out) <- NULL
  return(sdata_out)
}


# exportPEAK ------------------------------------------------------------------------------------------------------
#' @title Export centWave measures for a peak-of-interest in every sample in which it was detected.
#'
#' @param peakid \code{numeric} corresponding to peak ID in the final template.
#' @param peaks_vals_samples \code{list} with centWave measures for every sample, with length equal to sample number.
#'
#' @return Function return a \code{data.frame} with centWave measures for a peak-of-interest, for every sample in the provided peaks_vals_samples list.
#' 
exportPEAK <- function(peakid, peaks_vals_samples) {
  vals <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo")
  peak_all <- data.frame()
  ## for every sample
  for (n in 1:length(peaks_vals_samples)) {
    sdata <- peaks_vals_samples[[n]]
    peak_n <- sdata[match(peakid, sdata$peakid), vals]
    peak_all <- rbind(peak_all, peak_n, make.row.names = FALSE)
  }
  return(peak_all)
}  


# getPEAKmedians --------------------------------------------------------------------------------------------------
#' @title Get median centWave measures for a peak-of-interest across all samples in which it was detected.
#'
#' @param peak_n \code{data.frame} with centWave measures for a peak-of-interest in every sample
#' @param values \code{character} with centWave measures for which medians should be estimated.
#' Default set to: c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")
#'
#' @return Function return a \code{data.frame} with median centWave measures for a peak-of-interest across all samples.
#' 
getPEAKmedians <-
  function(peak_n,
           values = c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")) {
    medians <-
      apply(peak_n[, match(values, colnames(peak_n))], 2, median, na.rm = TRUE)
    npeaks <- length(which(!is.na(peak_n$mz)))
    peak_nmat <- as.data.frame(cbind(t(medians), npeaks))
    return(peak_nmat)
  }

