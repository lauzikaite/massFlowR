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
#' @param files \code{character} with paths to mzML file(s) to be processed.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param cwt \code{CentWaveParam} class object with parameters for centWave-based peak-picking.
#' @param ncores \code{numeric} defining number of cores to use for parallelisation. Default set to 1 for serial implementation.
#' @param thr \code{numeric} defining correlation coefficient threshold, above which peak pairs will be considered as correlated.
#'
#' @return For each LC-MS file, function writes a table with picked-peaks and their peak-groups into separate files, in the defined \code{out_dir} directory.
#' 
#' @export
#'
groupPEAKS <- function(files, out_dir, cwt, ncores = 1, thr = 0.95) {

  if (missing(files)) { 
    stop("'files' must be specified!")
  }
  if (missing(out_dir)) {
    stop("'out_dir' must be specified!")
  }
  if (missing(cwt)) {
    stop("'cwt' has to be specified!")
  } else {
    if (class(cwt) != "CentWaveParam") { 
      stop("'cwt' has to be 'CentWaveParam' object!") 
    }
    cwt@verboseColumns <- TRUE ## verboseColumns must be TRUE to output column "scpos"
  }
  if (!is.numeric(ncores)) {
    stop("'ncores' has to be numeric value!")
  } else {
    if (ncores < 1) {
      stop("'ncores' must be set to 1 (serial performance), or higher!")
    }
  }

  if (ncores > 1) {
    message("'ncores' set to ", ncores)
    nfiles <- length(files)
    ## get number of processes across which files will be divided
    nproc <- ceiling(nfiles/ncores)

    ## in each new process, start a new cluster with selected number of cores
    for (iproc in seq(0, (nproc -1), by = 1)) {
      
      ## extract files to be processed in this process
      files_proc_first <- 1 + (iproc*ncores) 
      files_proc_last <- ncores + (iproc*ncores)
      files_proc_last <- ifelse(files_proc_last <= nfiles, files_proc_last, nfiles) 
      files_proc <- files[files_proc_first:files_proc_last]
  
      ## run selected files across the cluster
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      foreach::foreach(f = files_proc, .inorder = TRUE) %dopar% groupPEAKS_paral(f = f,
                                                                               out_dir = out_dir,
                                                                               cwt = cwt,
                                                                               thr = thr)
      parallel::stopCluster(cl)
      message(files_proc_last, " out of ", nfiles, " files were processed.")
    }
    # ## use DoparParam for foreach implementation
    # BiocParallel::bplapply(files,
    #                        FUN = groupPEAKS_paral,
    #                        out_dir = out_dir,
    #                        cwt = cwt,
    #                        thr = thr,
    #                        BPPARAM = BiocParallel::DoparParam())
    ## close parallel backend
    # parallel::stopCluster(cl)
  } else {
    ## serial implementation
    message("'ncores' not defined. Running in serial mode ...")
    lapply(files,
           FUN = groupPEAKS_paral,
           out_dir = out_dir,
           cwt = cwt,
           thr = thr)
  }
  message("Peak-groups for all files were succesfully generated.")
}

# groupPEAKS_paral ------------------------------------------------------------------------------------------------------
groupPEAKS_paral <- function(f, out_dir, cwt, thr) {
  fname <- strsplit(basename(f), split = "[.]")[[1]][1]
  
  ## use try to catch mzML reading error that occurs on macOS
  raw <- NULL
  while (is.null(raw)) {
    raw <- try(MSnbase::readMSData(f, mode = "onDisk"),
               silent = TRUE)
    if (class(raw) == "try-error") {
      message("readMSData fail. failing mzML file: ", fname)
      message("reruning file ...")
      raw <- NULL
    }
  }
  pks <- pickPEAKS(raw = raw, cwt = cwt, fname = fname)
  eic <- extractEIC(raw = raw, pks = pks)
  groupPEAKSspec(pks = pks, eic = eic, out_dir = out_dir, fname = fname, thr = thr)
  }


