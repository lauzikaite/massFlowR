#' @title Get peak-groups in every LC-MS datafiles
#' @description Function performs peak-picking, based on package \emph{xcms} functionality.
#' Picked-peaks are then used to build groups of correlated, co-eluting peaks.
#'
#' @param files A \code{character} with paths to files to be processed.
#' @param out_dir A \code{character} specifying desired directory for output.
#' @param cwt A \code{CentWaveParam} class object with parameters for centWave-based peak-picking.
#' @param match A \code{numeric} defining the number of scans for co-eluting peaks extraction.
#' @param pearson A \code{logical} whether Pearson Correlation should be used. For \code{pearson = FALSE}, Spearman correlation method will be used.
#' @param thr A \code{numeric} defining correlation coefficient threshold, above which peak pairs will be considered as correlated.
#' @param plot A \code{logical}. For \code{plot = TRUE}, a network graph for each peak in the table will be saved as a png file in the out_dir directory.
#' @param clean A \code{logical} whether one-peak peak-groups should be removed (default is TRUE).
#' @param bpparam A \code{BiocParallel} parameter object to control how and if parallel processing should be performed. Such object can be created by the \emph{SerialParam}, \emph{MulticoreParam} or \emph{SnowParam} functions from the \emph{BiocParallel} package.
#'
#' @return For each LC-MS file, function writes a table with picked-peaks and their peak-groups into separate files, in the defined \code{out_dir} directory.
#' @export
#'
groupPEAKS <- function(files, out_dir, cwt, match = 1, pearson = TRUE, thr = 0.95, plot = FALSE, clean = TRUE, bpparam) {

  if (missing(files)) { stop("'files' must be specified!") }
  if (missing(out_dir)) { stop("'out_dir' must be specified!") }
  if (missing(cwt)) { stop("'cwt' has to be specified!") } else {
    if (class(cwt) != "CentWaveParam") { stop("'cwt' has to be 'CentWaveParam' object!") }
  }

  message("Apex matching window: ", match, " SCPOS")
  message("Correlation estimation: ", ifelse(pearson == TRUE, "Pearson", "Spearman"))

  ## if paral workers are not defined, use the default backend
  if (missing(bpparam)) { bpparam <-  BiocParallel::bpparam() }

  BiocParallel::bplapply(files,
                         FUN = groupPEAKS_paral,
                         out_dir = out_dir,
                         cwt = cwt,
                         match = match,
                         pearson = pearson,
                         thr = thr,
                         plot = plot,
                         clean = clean,
                         BPPARAM = bpparam)

  message("Peak-groups for all files were succesfully generated.")

}

####---- function for bplapply call, not to export
groupPEAKS_paral <- function(f, out_dir, cwt, match, pearson, thr, plot, clean) {

  fname <- strsplit(basename(f), split = "[.]")[[1]][1] # remove filename encoding for simplicity
  raw <- MSnbase::readMSData(f, mode = "onDisk")
  pks <- pickPEAKS(raw = raw, cwt = cwt, fname = fname, out_dir = out_dir)
  eic <- extractEIC(raw = raw, pks = pks)
  groupPEAKSspec(pks = pks, eic = eic, out_dir = out_dir, fname = fname, pearson = pearson, match = match, thr = thr, plot = plot, clean = clean)

}

