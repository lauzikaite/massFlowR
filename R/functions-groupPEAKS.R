# groupPEAKS ------------------------------------------------------------------------------------------------------
#' @aliases groupPEAKS
#' 
#' @title Get peak-groups representing different chemical spectra in each LC-MS datafile
#' 
#' @description Function performs peak-picking and peak-grouping for each LC-MS datafile independently.
#' 
#' @details Function performs peak-picking using the \emph{centWave} algorithm from package \emph{xcms}.
#' Picked-peaks are then grouped into chemical spectra. For each peak in the sample:
#' \itemize{
#' \item Co-eluting peaks are found.
#' \item EIC correlation between all co-eluting peaks is performed.
#' \item Network of peaks with high EIC correlation (with coefficients above the selected threshold) is built.
#' }
#' Co-eluting peaks within a network of high EIC correlation originate from the same chemical compound and therefore form a chemical spectrum.
#'
#' @param file \code{character} with path to metadata csv file, which must incluce columns 'filename', 'run_order' and 'raw_filepath'.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param cwt \code{CentWaveParam} class object with parameters for centWave-based peak-picking.
#' @param ncores \code{numeric} defining number of cores to use for parallelisation. Default set to 1 for serial implementation.
#' @param thr \code{numeric} defining correlation coefficient threshold, above which peak pairs will be considered as correlated. Default set to 0.95.
#'
#' @return For each LC-MS file, function writes a table with picked-peaks and their peak-groups into separate files, in the defined \code{out_dir} directory.
#' Inputed metadata file is updated to include a column 'proc_filepath', specifying full path to processed csv files peak-group tables.
#' 
#' @export
#'
groupPEAKS <- function(file = NULL,
                       out_dir = NULL,
                       cwt,
                       ncores = 2,
                       thr = 0.95) {
  if (is.null(file)) {
    stop("'file' is required")
  }
  if (!file.exists(file)) {
    stop("incorrect filepath for 'file' provided")
  }
  if (is.null(out_dir)) {
    stop("'out_dir' is required")
  }
  if (!dir.exists(out_dir)) {
    stop("incorrect filepath for 'out_dir' provided")
  }
  if (missing(cwt) | class(cwt) != "CentWaveParam") {
    stop("'cwt' has to be a 'CentWaveParam' object")
  } 
  if (ncores < 1 | !is.numeric(ncores)) {
    warning("'ncores' was not correctly set. Switching to ncores = 1 (serial performance)")
  }
  ## register paral backend
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  ## check provided metadata
  samples <- read.csv(file, header = T, stringsAsFactors = F)
  req_cnames <- c("filename",
                  "run_order",
                  "raw_filepath")
  if (any(!req_cnames %in% names(samples))) {
    stop("'files' table must contain columns: ", paste0(req_cnames, collapse = ", "))
  }
  ## verboseColumns must be TRUE to output column "scpos"
  cwt@verboseColumns <- TRUE
  ## mzdiff was >= 0 to avoid duplicate peaks in peak-groups
  if (cwt@mzdiff < 0) {
    message("mzdiff of ", cwt@mzdiff, " was selected. Switching mzdiff to 0 ..." )
    cwt@mzdiff <- 0
  }
  ## list raw files to be processed (only list files that have not been peak-grouped in the out_dir yet)
  samples$proc_filepath <- paste0(file.path(out_dir, samples$filename), "_peakgrs.csv")
  rfiles <- samples$raw_filepath[which(!file.exists(samples$proc_filepath))]

  while (length(rfiles) > 0) {
    
    ## create named vector to store processing result for every file
    nfiles <- length(rfiles)
    result <- setNames(vector('list', nfiles), nm = rfiles)
    
    if (ncores > 1) {
      
      ## get number of processes across which files will be divided
      nproc <- ceiling(nfiles/ncores)

      ## in each new process, start a new cluster with selected number of cores
      for (iproc in seq(0, (nproc -1), by = 1)) {
        
        ## extract files to be processed in this process
        files_proc_first <- 1 + (iproc*ncores) 
        files_proc_last <- ncores + (iproc*ncores)
        files_proc_last <- ifelse(files_proc_last <= nfiles, files_proc_last, nfiles) 
        files_proc <- rfiles[files_proc_first:files_proc_last]
    
        ## run selected files across the cluster
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
        result[files_proc_first:files_proc_last] <- foreach::foreach(
          f = files_proc,
          .inorder = TRUE,
          .errorhandling = "pass"
        ) %dopar% tryCatch({
          groupPEAKS_paral(
            f = f,
            out_dir = out_dir,
            cwt = cwt,
            thr = thr
          )
        },
        error = function(err) {
          return(list(fname = f, status = "FAILED", error = err$message))
        })
        parallel::stopCluster(cl)
        message(files_proc_last, " out of ", nfiles, " files were processed.")
      }
    } else {
      ## serial implementation
      result <- lapply(setNames(rfiles, rfiles),
                       function(f) {
                         tryCatch({
                           groupPEAKS_paral(
                             f = f,
                             out_dir = out_dir,
                             cwt = cwt,
                             thr = thr
                             )
                       },
                       error = function(err) {
                         return(list(
                           fname = f,
                           status = "FAILED",
                           error = err$message
                         ))
                       })})
    }
  
  ## process results
  # identify files that failed and update list of files
  result_status <- sapply(result, "[[", "status")
  rfiles <- names(result[which(result_status == "FAILED")])
  
  if (length(rfiles) > 0) {
    message(" \n", length(rfiles), " files failed. Processing failed files ...")
    message("Generated error messages for failed files:\n")
    message(paste0(sapply(result[rfiles], "[[", "error")))
  }
  }
  write.csv(samples, file = file, row.names = FALSE)
  message("Peak-groups for all files were succesfully generated.")
}

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
    groupPEAKSspec(pks = pks, eic = eic, out_dir = out_dir, fname = fname, thr = thr)
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
