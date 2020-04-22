# getRSE ------------------------------------------------------------------
#' @title Estimate RSE for a polynomial model with given peaks combination 
#'
#' @description Function build a polynomial model using the RT values for the target and its matching peaks combination. The residual standard error is returned.
#' 
#' @param combs \code{integer} with matching peaks peakids.
#' @param mat_dt \code{data.frame} with all peaks that belong to a given sub-peak.
#' @param t_peakid \code{integer} with the current target peak's peakid.
#'
#' @return Function returns RSE value for the model built with the given peaks combination.
#' 
getRSE <- function(combs, mat_dt, t_peakid) {
  combs_dt <- subset(mat_dt, peakid %in% combs | peakid == t_peakid)
  if (length(combs) > 1) {
    lmodel <- lm(rt ~ poly(sample, 2) + peakid, data = combs_dt)
  } else {
    lmodel <- lm(rt ~ poly(sample, 2), data = combs_dt)
  }
  return(round(summary(lmodel)$sigma, digits = 4))
}


# checkTARGETpeak ---------------------------------------------------------
#' @title Check a given target peak and all of its matching peaks.
#' 
#' @description Function evaluates all of the matching peaks for a given target peak.
#' Function performs polynomial model fitting using RT values. 
#' 
#' @param t_peakid \code{integer} with the current target peak's peakid.
#' @param target_dt \code{data.frame} with all peaks that belong to the target peak.
#' @param mat_all  \code{data.frame} with all matching peaks. 
#'
#' @return
checkTARGETpeak <- function(t_peakid, target_dt, mat_all) {
  t_peakid_dt <- subset(target_dt, peakid == t_peakid)
  
  ## only if target peak has atleast >10 samples (temporal threshold)
  if (nrow(t_peakid_dt) < 10) {
    return()
  }
  
  ## if target peak matches any peaks
  t_peakid_mat <- subset(mat_all, target_peakid == t_peakid)
  m_peakids <- unique(t_peakid_mat$peakid)
  if (length(m_peakids) == 0) {
    return()
  }
  t_peakid_mat_dt <- data.frame(
    sample = t_peakid_dt$"sample",
    rt = t_peakid_dt$"rt",
    # peakid = rep("target", times = nrow(t_peakid_dt)),
    peakid = rep(t_peakid, times = nrow(t_peakid_dt)),
    stringsAsFactors = FALSE
  )
  
  # compare completion patterns ---------------------------------------------
  for (m_peakid in m_peakids) {
    m <- subset(t_peakid_mat, peakid == m_peakid, select = c("sample", "rt", "peakid"))
    m_t <- unique(c(t_peakid_mat_dt$sample, m$sample))
    score <- sum(length(c(setdiff(t_peakid_mat_dt$sample, m$sample), setdiff(m$sample, t_peakid_mat_dt$sample))),
                 length(intersect(t_peakid_mat_dt$sample, m$sample))) / length(m_t)
    if (abs(1 - score) > 1/3) {
      # remove matching peak based on low completion pattern
      m_peakids <- m_peakids[-which(m_peakids == m_peakid)]
    }
  }
  if (length(m_peakids) == 0) {
    return()
  }
  
  # check for multimodality -------------------------------------------------
  ## if target peak itself have distinct RT populations
  dE <- diptest::dip.test(t_peakid_mat_dt$rt, simulate.p.value = TRUE, B = 1000)
  ## find local maxima in the RTs of the target peak
  dn <- density(t_peakid_mat_dt$rt)
  dn_max <- dn$x[which(diff(sign(diff(dn$y)))==-2)+1]
  ## if multiple maxima points found, the bimodial mode is no good
  if (dE$p.value <= 0.01 & length(dn_max) == 2) {
    if (diff(dn_max) > rt_err) {
      ## find standard deviation around the maxima points
      t_peakid_mat_dt$sub_peak <- sapply(t_peakid_mat_dt$rt, function(rt) {
        which(abs(dn_max - rt) == min(abs(dn_max - rt)))
      })
      dn_max_err <- lapply(seq(length(dn_max)), function(s_peak) {
        rts <- subset(t_peakid_mat_dt, sub_peak == s_peak)$rt
        std <- abs(sd(rts)) * 3
        c(dn_max[s_peak] - std, dn_max[s_peak] + std)
      })
      ## assign matching peaks to the closest sub-peak
      for (m_peakid in m_peakids) {
        m <- subset(t_peakid_mat, peakid == m_peakid, select = c("sample", "rt", "peakid"))
        m$sub_peak <- sapply(m$rt, function(rt) {
          ## to which sub-peak this RT is closest?
          subpk <- which(abs(dn_max - rt) == min(abs(dn_max - rt)))
          ## is the distance to this sub-peak still within the error range?
          if (all(rt >= dn_max_err[[subpk]][[1]] & rt <= dn_max_err[[subpk]][[2]])) {
            subpk
          } else {
            NA
          }
        })
        ## do not add the peak if it is not assigned to either of the sub-peaks in none of the samples
        if (all(is.na(m$sub_peak))) {
          # remove matching peak
          m_peakids <- m_peakids[-which(m_peakids == m_peakid)]
        } else {
          t_peakid_mat_dt <- rbind(t_peakid_mat_dt, m)
        }
      }
      if (length(m_peakids) == 0) {
        return()
      }
    }
  } else {
    t_peakid_mat_dt <- rbind(t_peakid_mat_dt,
                             subset(t_peakid_mat, peakid %in% m_peakids)[ , c("sample", "rt", "peakid")])
    t_peakid_mat_dt$sub_peak <- 1
  }
  
  # model drifts ------------------------------------------------------------
  ## compare each sub-peak with its assigned matches
  sub_peaks <- na.omit(unique(t_peakid_mat_dt$sub_peak))
  t_peakid_mat_dt$joined_sub_peak <- NA
  for (s_peak in sub_peaks) {
    s_peak_dt <- subset(t_peakid_mat_dt, sub_peak == s_peak)
    s_peak_m_peakids <- unique(subset(s_peak_dt, peakid != t_peakid)$peakid)
    if (length(s_peak_m_peakids) > 0) {
      ## regression model for just the target peak
      s_peak_rse <- round(summary(lm(rt ~ poly(sample, 2),
                                     data = subset(s_peak_dt, peakid == t_peakid)))$sigma, digits = 4)
      ## regression models with incrementally added combinations of matching peaks
      for (m in seq(from = length(s_peak_m_peakids), to = 1)) {
        if (length(s_peak_m_peakids) > 1) {
          combination <- gRbase::combnPrim(x = s_peak_m_peakids, m = m, simplify = FALSE)
        } else {
          combination <- list(s_peak_m_peakids)
        }
        combination_rse <- lapply(combination, FUN = getRSE, mat_dt = s_peak_dt, t_peakid = t_peakid)
        if (any(combination_rse <= s_peak_rse)) {
          t_peakid_mat_dt$joined_sub_peak[which(
            t_peakid_mat_dt$peakid %in% c(t_peakid, combination[[which.min(combination_rse)]]) &
              t_peakid_mat_dt$sub_peak == s_peak)] <- s_peak
          break()
        }
      }
    } else {
      t_peakid_mat_dt$joined_sub_peak[which(t_peakid_mat_dt$sub_peak == s_peak)] <- s_peak
    }
  }
  
  ## return full information on how target and matches were split into sub-peaks
  t_peakid_mat_dt$joined_peakid <- t_peakid
  t_peakid_mat_dt$joined_peakid[which(is.na(t_peakid_mat_dt$joined_sub_peak))] <- NA
  return(t_peakid_mat_dt)
}