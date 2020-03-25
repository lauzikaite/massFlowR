# modelPEAKS -------------------------------------------------------------------------------------------------------
#' @title Function models intensity integration values for missing peaks.
#'
#' @description Function models intensity integration values for each of the missing peaks and for each of the samples in which it was missed.
#' Modelling relies on a smooth spline fit to the values of the corresponding detected peaks.
#'
#' @param p \code{numeric} indicating peak IDs.
#' @param vars \code{character} with centWave values to be modelled, e.g. c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax").
#' @param object \code{massFlowTemplate} class object.
#'
#' @return Function returns a list of data.frames - a data.frame for each peak supplied with the argument p.
#'
modelPEAKS <- function(p, vars, object) {
  peak <- object@peaks[[match(p, names(object@peaks))]]
  pred_vals <- lapply(vars,
    FUN = smoothVALUE,
    x = 1:nrow(peak),
    peak = peak
  )
  peak_vals <- setNames(as.data.frame(do.call("cbind", pred_vals)),
    nm = vars
  )
  return(peak_vals)
}

# smoothVALUE -------------------------------------------------------------------------------------------------------
#' @title Helper function to model intensity integration values for missing peaks.
#'
#' @description Helper function is used by the \code{\link{modelPEAKS}} function to model intensity integration values for missing peaks.
#' Function fit a smoothing spline to centWave values reported for the peak-of-interest in all samples in the study.
#' For samples in which peak was not detected, centWave values are predicted.
#'
#' @param var \code{character} with a single centWave value to be modelled, e.g. "mz"
#' @param x \code{numeric} representing indeces of peaks for the supplied 'peak' data.frame.
#' @param peak \code{data.frame} with all centWave values (columns) for the peak-of-interest for all samples (rows). For samples in which peak as not detected, NAs must be listed.
#' @param lambda \code{numeric} representing smoothing parameter lambda. Set to 10 by default, which corresponds to number of study samoples in between two quality control samples in standard experimental protocols.
#'
#' @return Function returns centWave values for the peak-of-interest, which includes predicted values for samples in which it was not originally detected.
#'
smoothVALUE <- function(var, x, peak, lambda = 10) {
  y <- peak[, match(var, names(peak))]
  present <- which(!is.na(y))
  if (length(x[present]) < 4 ) {
    ## spline will fail with an error message otherwise
    return(y)
  }
  mod <- stats::smooth.spline(
    x = x[present],
    y = y[present],
    lambda = lambda,
    cv = TRUE
  )
  predicted <- stats::predict(mod, x = x[-present])$y
  y[is.na(y)] <- predicted
  return(y)
}

# extractMISS  -------------------------------------------------------------------------------------------------------
#' @title Extract samples with missed peaks, for which intensity have to be re-integrated.
#'
#' @description Function builds a list with mz and rt regions for intensity re-integration.
#' This is a data-reduction step, helping to reduce the size of the variables being imputed to \code{\link{fillSAMPLE}} function for intensity re-integration.
#'
#' @param s \code{numeric} corresponding to sample number in the study.
#' @param values \code{list} with centWave values for each sample.
#' @param peaks_mod \code{list} with real and modelled centWave values for each peak.
#' @param valid \code{data.frame} with validated aligned peaks.
#'
#' @return Function returns a list with a data.frame for each sample in the study.
#' Each data.frame contains predicted mz and rt regions only for those peaks that were missed in the corresponding sample.
#'
extractMISS <- function(s, values, peaks_mod, valid) {
  ## for each missing peak, take the modelled mz/rt values and return in a new df
  sample <- values[[s]]
  peaks_miss <- which(is.na(sample$mz))
  if (length(peaks_miss) == 0) {
    return(data.frame())
  }
  peaks_mod_miss <- lapply(peaks_mod[peaks_miss], function(p) {
    p[s, ]
  })
  peaks_mod_miss <- do.call("rbindCLEAN", peaks_mod_miss)
  peaks_mod_miss$peakid <- valid$peakid[peaks_miss]
  return(peaks_mod_miss)
}

# fillSAMPLE -------------------------------------------------------------------------------------------------
#' @title Fill intensity values for peaks missed in the sample-of-interest.
#'
#' @param s \code{numeric} corresponding to the number of the sample-of-interest in the study.
#' @param sname \code{character} with full path to the raw LC-MS file for the sample-of-interest.
#' @param sdata \code{data.frame} with all centWave values for the sample-of-interest.
#' @param values \code{data.frame} with predicted mz and rt regions only for those peaks that were missed in the corresponding sample.
#'
#' @return Function returns updated \code{data.frame} with the centWave values for the sample-of-interest.
#'
fillSAMPLE <- function(s, sname, sdata, values) {
  if (validFILE(sname) == TRUE) {
    ## if missing peaks are present for this sample
    if (nrow(values) > 0) {
      raw <- readDATA(sname)
      ## extract intensity of scans corresponding to mz&rt regions for every peak
      peaks_int <- extractINT(raw,
        mz = values[, c("mzmin", "mzmax")],
        rt = values[, c("rtmin", "rtmax")]
      )
  
      ## get into and maxo values for each peak
      values_filled <- as.data.frame(setNames(replicate(n = ncol(values) + 2, numeric(), simplify = FALSE), nm = c(names(values), "maxo", "into")))
  
      for (p in 1:length(peaks_int)) {
        peak <- peaks_int[[p]]
        values_filled[p, names(values)] <- values[p, ]
  
        ## if scans were found for this peak
        if (nrow(peak) > 0) {
          ## get rt width
          rts <- unique(peak$rt)
          rt_range <- range(peak$rt)
          rt_width <- (rt_range[2] - rt_range[1]) / length(rts)
  
          ## get max int for each scan
          rt_max <- sapply(rts, function(rt) {
            max(peak[peak$rt == rt, "int"])
          })
          maxo <- max(rt_max)
  
          ## get into integration
          mz_max <- sapply(unique(peak$mz), function(mz) {
            max(peak[peak$mz == mz, "int"])
          })
          into <- sum(mz_max) * rt_width
          values_filled[p, c("maxo", "into")] <- c(maxo, into)
        } else {
          ## return default values - 0s
          values_filled[p, c("maxo", "into")] <- 0
        }
      }
      ## bind with centWave results
      ## retain original peakid order in sdata
      sdata[match(values_filled$peakid, sdata$peakid), colnames(values_filled)] <- values_filled
  }
    ## cleanup garbage
    gc(verbose = FALSE)
    return(sdata)
  } else {
    return(list(nsample = s, status = "FAILED", error = validFILE(sname)))
  }
}

# extractINT --------------------------------------------------------------
#' @title Extract raw signal values for defined mz and rt region.
#'
#' @description Helper function used by \code{\link{fillSAMPLE}} function to extract spectrum intensity values for a defined mz and rt window.
#'
#' @param raw \code{OnDiskMSnExp} class object.
#' @param rt \code{data.frame} with rtmin and rtmax values (columns) for peaks (rows) that have to be re-integrated in the sample-of-interest.
#' @param mz \code{data.frame} with mzmin and mzmax values (columns) for peaks (rows) that have to be re-integrated in the sample-of-interest.
#'
#' @return Function returns a \code{list} of data.frames with extracted raw LC-MS signal values.
#'
extractINT <- function(raw, rt, mz) {

  ## which scans are required to be filled
  raw_rt <- MSnbase::rtime(raw)
  scans_to_fill <- sort(unique(unlist(apply(rt, 1, function(r) {
    which(raw_rt >= r[1] &
      raw_rt <= r[2])
  }))))
  # if no scans
  if (length(scans_to_fill) == 0) {
    return(data.frame())
  }
  ## extract spectra and rt values for each scan
  raw_rt_to_fill <- raw[scans_to_fill]
  raw_spec <- MSnbase::spectra(raw_rt_to_fill)
  raw_spec_rt <- MSnbase::rtime(raw_rt_to_fill)

  ## filter by mz for each peak in the list
  peaks_int <- vector("list", nrow(rt))
  for (p in 1:nrow(rt)) {
    ## extract scans for the single peak
    peak_scans <- which(raw_spec_rt >= rt[p, 1] &
      raw_spec_rt <= rt[p, 2])
    if (length(peak_scans) == 0) {
      peak_int <- data.frame()
    } else {
      ## extract corresponding spectra values
      peak_scans <- raw_spec[peak_scans]
      peak_mz <- mz[p, ]
      peak_int <- lapply(peak_scans, function(scan) {
        suppressWarnings(
          scan_mz <- MSnbase::filterMz(scan, peak_mz)
        )
        if (!scan_mz@peaksCount) {
          NULL
        } else {
          data.frame(rt = scan_mz@rt, mz = scan_mz@mz, int = scan_mz@intensity)
        }
      })
      peak_int <- Filter(length, peak_int)
      if (length(peak_int) > 0) {
        peak_int <- do.call("rbindCLEAN", peak_int)
      } else {
        peak_int <- data.frame()
      }
    }
    peaks_int[[p]] <- peak_int
  }
  return(peaks_int)
}
