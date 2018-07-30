#' @title Obtain components for all datafiles
#'
#' @param files
#' @param out_dir
#' @param cwt
#' @param match
#' @param pearson
#' @param thr
#' @param plot
#' @param clean
#' @param workers
#'
#' @return
#' @export
#'
#' @examples
getCOMPS <- function(files, out_dir, cwt, match = 1, pearson = TRUE, thr = 0.95, plot = FALSE, clean = TRUE, bpparam) {

  if (missing(out_dir)) { stop("'out_dir' must be specified!") }
  if (missing(cwt)) { stop("'cwt' has to be specified!") }
  if (class(cwt) != "CentWaveParam") { stop("'cwt' has to be 'CentWaveParam' object!") }

  message("Apex matching window: ", match, " SCPOS")
  message("Correlation estimation: ", ifelse(pearson == TRUE, "Pearson", "Spearman"))

  ## if paral workers are not defined, use the default backend
  if (missing(bpparam)) { bpparam <-  BiocParallel::bpparam() }

  BiocParallel::bplapply(X = files,
                         FUN = getCOMPS_paral,
                         out_dir = out_dir,
                         cwt = cwt,
                         match = match,
                         pearson = pearson,
                         thr = thr,
                         plot = plot,
                         clean = clean,
                         BPPARAM = bpparam)

  message("Components were succesfully generated.")

}

####---- function for bplapply call, not to export
getCOMPS_paral <- function(f, out_dir, cwt, match, pearson, thr, plot, clean) {

  fname <- strsplit(basename(f), split = "[.]")[[1]][1] # remove filename encoding for simplicity
  raw <- MSnbase::readMSData(f, mode = "onDisk")
  pks <- pickPEAKS(raw = raw, cwt = cwt, fname = fname, out_dir = out_dir)
  eic <- extractEIC(raw = raw, pks = pks)
  buildCOMPS(pks = pks, eic = eic, out_dir = out_dir, fname = fname, pearson = pearson, match = match, thr = thr, plot = plot, clean = clean)

}

