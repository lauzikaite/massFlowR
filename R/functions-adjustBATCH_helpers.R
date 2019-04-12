checkCOMPOUND <- function(rc, val_pks, roi_end, roi_start, data_end, data_start) {
  ## find corresponding peak in the end sample
  peak_end <- findMRpeak(data = data_end, roi =  roi_end[rc, ])
  if (nrow(peak_end) > 0) {
    ## check if in the final template
    if (peak_end$tmp_peakid %in% val_pks$peakid) {
      ## find corresponding peak in the start sample
      peak_start <- findMRpeak(data = data_start, roi =  roi_start[rc, ])
      if (nrow(peak_start) > 0) {
        ## take peaks mz&rt values and get differences
        return(peak_start$rt - peak_end$rt)
      }
    }
  }
  return(NULL)
}

  
  
 
findMRpeak <- function(data, roi) {
  peak <- data[which(
    data$mzmin >= roi$mzMin &
      data$mzmax <= roi$mzMax &
      data$rtmin >= roi$rtMin &
      data$rtmax <= roi$rtMax
  ),]
  if (nrow(peak) > 1) {
    close_by_mz <- which.min(peak$mz - roi$mz)
    peak <- peak[close_by_mz, ]
  }
  return(peak)
}
 