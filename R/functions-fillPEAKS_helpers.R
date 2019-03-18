# modelPEAKS -------------------------------------------------------------------------------------------------------
modelPEAKS <- function(p, vars, object) {
  peak <- object@peaks[[match(p, names(object@peaks))]]
  pred_vals <- lapply(vars,
                      FUN = smoothVALUE,
                      x = 1:nrow(peak),
                      peak = peak)
  peak_vals <- setNames(as.data.frame(do.call("cbind", pred_vals)),
                        nm = vars)
  # peak_vals$sample_order <- 1:nrow(peak)
  return(peak_vals)
}

# smoothVALUE -------------------------------------------------------------------------------------------------------
smoothVALUE <- function(var, x, peak, lambda = 10) {
  y <- peak[,match(var, names(peak))]
  present <- which(!is.na(y))
  mod <- smooth.spline(x = x[present],
                       y = y[present],
                       lambda = lambda,
                       cv = TRUE)
  predicted <- predict(mod, x = x[-present])$y
  y[is.na(y)] <- predicted
  return(y)
}

# extractMISS  -------------------------------------------------------------------------------------------------------
extractMISS <- function(s, values, peaks_mod, valid) {
  ## for each missing peak, take the modelled mz/rt values and return in a new df
  sample <- values[[s]]
  peaks_miss <- which(is.na(sample$mz))
  peaks_mod_miss <- lapply(peaks_mod[peaks_miss], function(p) {
    p[s,]
  })
  peaks_mod_miss <- do.call("rbindCLEAN", peaks_mod_miss)
  peaks_mod_miss$peakid <- valid$peakid[peaks_miss]
  return(peaks_mod_miss)
}

# fillSAMPLE -------------------------------------------------------------------------------------------------
fillSAMPLE <- function(s, sname, sdata, values) {
  
  if (validFILE(sname) == TRUE) {
    raw <- readDATA(sname)
    ## extract intensity of scans corresponding to mz&rt regions for every peak
    peaks_int <- extractINT(raw,
                            mz = values[, c("mzmin", "mzmax")],
                            rt = values[, c("rtmin", "rtmax")])
    ## sanity check for now
    if (length(peaks_int) != nrow(values)) {
      stop("fillSAMPLE fails")
    }
    
    ## get into and maxo values for each peak
    values_filled <- data.frame()
    for (p in 1:length(peaks_int)) {
      peak  <- peaks_int[[p]]
      values_filled[p, names(values)] <- values[p,]
      
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
        ## return default values
        values_filled[p, c("maxo", "into")] <- 0
      }
    }
    
    ## bind with centWave results
    ## retain original peakid order in sdata
    sdata[match(values_filled$peakid, sdata$peakid), colnames(values_filled)] <- values_filled
  
    ## cleanup garbage
    gc(verbose = FALSE)
    return(sdata)
  } else {
    return(list(nsample = s, status = "FAILED", error = validFILE(sname)))
  }
} 


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
        if(!scan_mz@peaksCount) {
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

