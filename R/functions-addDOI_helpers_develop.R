matchPEAK2 <- function(peakgr, doi, tmp) {
  mat <- lapply(peakgr,
                FUN = matchPEAK3,
                doi = doi,
                tmp = tmp)
  if (length(mat) != length(peakgr)) {
    stop("matchPEAK3 returns mistake!")
  }
  return(mat)
  
}

matchPEAK3 <- function(peak, doi, tmp) {
  target_peak <- doi[peak, c("mz_l", "mz_h", "rt_l", "rt_h")]
  matches_mz <-
    tmp[which(
      dplyr::between(tmp$mz_l, target_peak$mz_l, target_peak$mz_h) |
        dplyr::between(tmp$mz_h,  target_peak$mz_l, target_peak$mz_h)
    ), ]
  matches_mz_rt <-
    which(
      dplyr::between(matches_mz$rt_l, target_peak$rt_l, target_peak$rt_h) |
        dplyr::between(matches_mz$rt_h, target_peak$rt_l, target_peak$rt_h)
    )
  if (length(matches_mz_rt) > 0) {
    tmp_ind <- as.numeric(rownames(matches_mz)[matches_mz_rt])
  } else {
    tmp_ind <- 0
  }
  return(tmp_ind)
}

comparePEAKGR <- function(peakgr, target, doi, mat, tmp, bins) {
  ## extract peaks for this target peakgr
  target_peaks <- doi[target[[peakgr]], c("mz", "into")]
  
  ## extract matched peaks
  ## 1. extract all tmp peaks, that belong to the same clusters as the peaks matched directly by mz/rt
  matched_ind <- unlist(mat[[peakgr]])
  
  matched_peaks <- tmp[matched_ind, c("peakgr", "peakgrcls")]
  matched_peakgrcls <- unique(matched_peaks$peakgrcls)
  matched_peakgrs <- unique(matched_peaks$peakgr)
  matched_peaks_all <-
    tmp[which(tmp$peakgrcls %in% matched_peakgrcls), c("mz", "into", "peakgr", "peakgrcls")]
  
  ## get total mz range between the target peaks and all matched template peaks
  ## generate mz bins using all peaks
  all_peaks <- c(target_peaks$mz, matched_peaks_all$mz)
  breaks <-
    data.frame(breaks = seq(
      from = min(all_peaks),
      to = max(all_peaks),
      by = bins
    ))
  breaks$bin <- 1:nrow(breaks)
  
  ## build target vector of the target peak-group spectrum
  target_peaks$bin <-
    findInterval(x = target_peaks$mz, breaks$breaks)
  ## remove duplicating duplicating bins
  target_vec <-
    target_peaks[which(!duplicated(target_peaks$bin)), c("bin", "into")]
  target_vec <- dplyr::full_join(breaks,
                                 target_vec,
                                 by = c("bin"))
  target_vec <- ifelse(is.na(target_vec$into), 0, target_vec$into)
  target_vec <- scaleVEC(target_vec)
  
  ## for every matched peakgr, find cosine
  matched_peakgrs_cos <- sapply(
    matched_peakgrs,
    FUN = getCOS2,
    target_vec = target_vec,
    breaks = breaks,
    matched_peaks_all = matched_peaks_all
  )
  matched_peakgrs_cos <- data.frame(
    target_peakgr = names(target)[[peakgr]],
    matched_peakgr = matched_peakgrs,
    cos = matched_peakgrs_cos
  )
  
  return(matched_peakgrs_cos)
}

getCOS2 <- function(peakgr,
                    target_vec,
                    breaks,
                    matched_peaks_all) {
  ## build vector from matched peaks
  matched_vec <-
    matched_peaks_all[which(matched_peaks_all$peakgr == peakgr), ]
  matched_vec$bin <- findInterval(x = matched_vec$mz, breaks$breaks)
  ## remove duplicating duplicating bins
  matched_vec <-
    matched_vec[which(!duplicated(matched_vec$bin)), c("bin", "into")]
  matched_vec <- dplyr::full_join(breaks,
                                  matched_vec,
                                  by = c("bin"))
  matched_vec <-
    ifelse(is.na(matched_vec$into), 0, matched_vec$into)
  matched_vec <- scaleVEC(matched_vec)
  
  ## find the cosine of the angle
  cos_angle <-
    (sum(target_vec * matched_vec))  / ((sqrt(sum(
      target_vec * target_vec
    )))  * (sqrt(sum(
      matched_vec * matched_vec
    ))))
  
  return(cos_angle)
}
