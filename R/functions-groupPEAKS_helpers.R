# groupPEAKS_paral ------------------------------------------------------------------------------------------------------
#' @title Perform peak-picking and grouping into Extracted Chemical Spectra for a single file
#' 
#' @description Function enabes parallel or serial implementation of peak-picking and grouping.
#' This function is called from inside \code{\link{groupPEAKS}}, for each sample in the experiment separately.
#' Function is a wrapper for multiple individual raw file processing steps, which depend on \code{MSnbase} functionality.
#' To overcome raw file reading into memory issues associated with MacOS, each function is within a while and try loop.
#' 
#' @param f \code{character} specifying absolute path to a single mzML/CDF file.
#' @param out_dir \code{character} specifying absolute path to directory for output.
#' @param cwt \code{CentWaveParam} class object with parameters for peak-picking. Object can be created by the \emph{xcms::CentWaveParam} function.
#' @param thr \code{numeric} defining Pearson correlation coefficient threshold, above which peak pairs will be considered as correlated. 
#'
#' @return Function returns a list with filename, status and error message (deafult is NULL, gets overwritten by tryCatch is function fails).
#' 
groupPEAKS_paral <- function(f, out_dir, cwt, thr) {
  fname <- strsplit(basename(f), split = "[.]")[[1]][1]
  
  if (validFILE(f) == TRUE) {
    raw <- readDATA(f = f)
    pks <- pickPEAKS(raw = raw, cwt = cwt, fname = fname)
    eic <- extractEIC(raw = raw, pks = pks)
    do_groupPEAKS(pks = pks, eic = eic, out_dir = out_dir, fname = fname, thr = thr)
    result <- list(fname = fname, status = "DONE", error = NULL)
  } else {
    result <- list(fname = fname, status = "FAILED", error = validFILE(f))
  }
  gc(verbose = FALSE)
  return(result)
}


# pickPEAKS -------------------------------------------------------------------------------------------------------
#' @aliases pickPEAKS
#'
#' @title Peak detection using the centWave method
#'
#' @description High resolution, centroid LC/MS datafiles are peak-picked using centWave algorithm, implemented by \emph{xcms} package [Tautenhahn 2008].
#' \code{pickPEAKS} is applied to a single datafile.
#'
#' @param cwt \code{CentWaveParam} class object with parameters for peak-picking. Object can be created by the \emph{xcms::CentWaveParam} function.
#' @param raw \code{OnDiskMSnExp} class object for the single LC-MS spectrum of interest. Such object can be created by the \emph{MSnbase::readMSData} function.
#' @param fname \code{character} specifying datafile name.
#'
#' @return Function returns a \code{DataFrame} object with picked-peaks and centWave details.
#'
pickPEAKS <- function(raw, cwt, fname) {
  
  message("Peak-picking file: ", fname, " ...")
  
  ## use try to catch mzML reading error that occurs on macOS
  res <- NULL
  while (is.null(res)) {
    res <- try(xcms::findChromPeaks(object = raw,
                                    param = cwt),
               silent = TRUE)
    if (class(res) == "try-error") {
      message("reruning file ...")
      res <- NULL
    }
  }
  pks <- data.frame(xcms::chromPeaks(res))
  
  ## order by intensity
  pks <- pks[order(pks$into, decreasing = T),]
  
  ## remove artefactural, duplicating peaks
  pks_unique <- unique(pks[,c("mz", "rt")])
  pks_clean <- lapply(1:nrow(pks_unique),
                      FUN = cleanPEAKS,
                      dt_unique = pks_unique,
                      dt = pks)
  pks_clean <- do.call("rbindCLEAN", pks_clean)
  
  ## assign each peak with id, based on its intensity order
  pks_clean$peakid <- 1:nrow(pks_clean)
  return(pks_clean)
}


# extractEIC ------------------------------------------------------------------------------------------------------
#' @title Obtain extracted ion chromatograms for all picked peaks
#'
#' @description Obtain extracted ion chromatograms (EIC) for all picked peaks, provided by \code{pickPEAKS} function.
#'
#' @param raw \code{OnDiskMSnExp} class object for the single datafile of interest.
#' @param pks \code{DataFrame} object generated with function \code{pickPEAKS}.
#'
#' @return Function returns a \code{list} with an EIC for eack peak in the \code{pks} table.
#' Individual peaks' EICs can be plotted using base graphics.
#'
extractEIC <- function(raw, pks) {
  eic <- xcms::chromatogram(
    raw,
    rt = data.frame(rt_lower = pks$rtmin,
                    rt_upper = pks$rtmax),
    mz = data.frame(mz_lower = pks$mzmin,
                    mz_upper = pks$mzmax)
  )
  clean_eic <- lapply(1:nrow(eic), function(ch) {
    MSnbase::clean(eic[ch,], na.rm = T)
  })
  return(clean_eic)
}

# do_groupPEAKS ------------------------------------------------------------------------------------------------------
#' @title Group peaks of a single LC-MS spectrum using EIC correlation
#'
#' @description Function builds peak-groups of co-eluting peaks that are correlated to each other in a single LC-MS spectrum.
#'
#' @param pks \code{DataFrame} containing the peak table, created by \emph{pickPEAKS} function.
#' @param eic \code{list} containing extracted ion chromatograms for each peak in the \code{pks} table.
#' @param out_dir \code{character} specifying directory where output data will be saved.
#' @param fname \code{character} specifying LC-MS filename.
#' @param thr \code{numeric} defining Pearson correlation coefficient threshold, above which peak pairs will be considered as correlated.
#' @param return \code{logical} whether to return the obtained peak-group table. Set to TRUE if working with function interactivily.
#'
#' @return Function returns an updated peak table.
#' Column \code{"peakgr"} contains peak-group ID for each peak.
#' Peaks not correlated to any of its co-eluting peaks are removed from the table.
#' Function also writes updated peak table in the specified directory.
#'
do_groupPEAKS <-
  function(pks, eic, out_dir, fname, thr, return = FALSE) {
    message("Building peak-groups for file: ", fname, "...")
    
    ## make a table to store peak-groups during grouping
    peakgroups <- pks[, c("peakid", "mz", "scpos", "into")]
    peakgroups$peakgr <- NA
    peakids <- peakgroups$peakid
    
    ## select POI (peak of interest) from peak table
    for (p in peakids) {
      ## if POI is already assigned to a group, skip to next
      if (!is.na(peakgroups[p, "peakgr"])) {
        next
      }
      ## find co-eluting peaks not asigned to a peakgronent yet
      ind <- c(peakgroups[p, "scpos"] - 1, peakgroups[p, "scpos"] + 1)
      co <- peakgroups[which(is.na(peakgroups$peakgr) &
                               (peakgroups$scpos >= ind[1] &
                                  peakgroups$scpos <= ind[2])),]
      co_ind <- co$peakid
      
      ## if POI doesn't have co-eluting peaks, skip to next
      if (length(co_ind) == 1 | !p %in% co_ind) {
        next
      }
      
      ## correlate the EIC of the peaks in the correlation matrix (cormat)
      cormat <- getCORmat(ind = co_ind)
      cormat$weight <- apply(cormat, 1, FUN = corEIC, eic = eic)
      
      ## if there isn't a single pair with cor above threshold, skip to next
      if (all(cormat$weight < thr)) {
        next
      }
      
      ## build network between co-eluting peaks and obtain communities
      coms <-
        buildGRAPH(pkg_cor = cormat,
                   cor_thr = thr,
                   plot = FALSE)
      
      ## extract community which includes the POI
      main_com <- coms[which(names(coms) == p)]
      
      ## extract all co-eluting peaks that are part of this community
      co_ind_com <- as.numeric(names(coms[which(coms == main_com)]))
      
      ## if main community only has the single peak, skip to next
      if (length(co_ind_com) == 1) {
        next
      }
      
      ## update table with assigned peakgroup IDs
      ## start with 1 if this is the first peakgroup to be assigned
      peakgr <-
        ifelse(all(is.na(peakgroups$peakgr)), 1, max(peakgroups$peakgr, na.rm = T) + 1)
      peakgroups[co_ind_com, "peakgr"] <- peakgr
    }
    
    ## removing peaks not assigned to any peak-group
    peakgroups <- peakgroups[which(!is.na(peakgroups$peakgr)), ]
    peakgroups <- merge(pks,
                        peakgroups[, c("peakid", "peakgr")],
                        by = c("peakid"), all = F)
    
    ## update peakids
    peakgroups$peakid <- 1:nrow(peakgroups)
    write.csv(
      peakgroups,
      file = paste0(out_dir, "/", fname, "_peakgrs.csv"),
      quote = F,
      row.names = F
    )
    
    if (return == TRUE) {
      return(peakgroups)
    } else {
      message(max(na.omit(peakgroups$peakgr)), " peak-groups built.")
    }
  }

# corEIC ------------------------------------------------------------------------------------------------------
#' @title Obtain EIC correlation between two peaks
#'
#' @description Function performs extracted ion chromatogram (EIC) correlation of two peaks.
#'
#' @param pair \code{matrix} with columns 'from' and 'to', specifying peak identifiers, which correspond to the indeces of a list of EICs.
#' @param eic \code{list} containing EIC for each picked-peak in a sample.
#'
#' @return Function returns EIC correlation coefficient.
#'
corEIC <- function(pair, eic) {
  x <- pair["from"]
  y <- pair["to"]
  rx <- MSnbase::rtime(eic[[x]])
  ry <- MSnbase::rtime(eic[[y]])
  common_scan <- base::intersect(rx, ry)
  if (length(common_scan) > 3) {
    ix <-
      as.numeric(MSnbase::intensity(eic[[x]])[which(rx %in% common_scan)])
    iy <-
      as.numeric(MSnbase::intensity(eic[[y]])[which(ry %in% common_scan)])
    cc <- cor(ix, iy, method = "pearson", use = "pairwise.complete.obs")
  } else {
    cc <- 0
  }
  ## negative coeficients would break graph generation
  cc <- ifelse(cc < 0, 0, cc)
  return(cc)
}