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
#' @param write
#' @param workers
#'
#' @return
#' @export
#'
#' @examples
getCOMPS <- function(files, out_dir, cwt, match = 1, pearson = TRUE, thr = 0.95, plot = FALSE, clean = TRUE, write = TRUE, workers = 2){

  message("Apex matching window: ", match, " SCPOS")
  message("Correlation estimation: ", ifelse(pearson == TRUE, "Pearson", "Spearman"))

  BPPARAM = BiocParallel::MulticoreParam(workers = workers, log = F, progressbar = T)

  BiocParallel::bplapply(X = files,
                         FUN = getCOMPS_paral,
                         out_dir = out_dir,
                         cwt = cwt,
                         match = match,
                         pearson = pearson,
                         thr = thr,
                         plot = plot,
                         clean = clean,
                         write = write,
                         BPPARAM = BPPARAM)

  message("Components were succesfully generated")

}

####---- function for bplapply call, not to export ----
getCOMPS_paral <- function(f, out_dir, cwt, match, pearson, thr, plot, clean, write) {

  if (missing(cwt)) { stop("'cwt' has to be specified!") }
  if (class(cwt) != "CentWaveParam") { stop("'cwt' has to be 'CentWaveParam' object!") }
  if(write == TRUE) { if(missing(out_dir)) { stop("'out_dir' must be specified!") } }

  fname <- gsub(".mzML", "", basename(f))
  raw <- MSnbase::readMSData(f, mode = "onDisk")
  pks <- pickPEAKS(raw = raw, cwt = cwt, fname = fname, out_dir = out_dir, write = write)

  eic <- extractEIC(raw = raw, pks = pks)

  # save(list = ls(envir = environment()), # must specify the environment, otherwise objects in the R_GlobalEnv will be saved
  #      file =  paste0(out_dir, fname, "_pks-eic.RData"))

  buildCOMPS(pks = pks, eic = eic, out_dir = out_dir, fname = fname, pearson = pearson, match = match, thr = thr, plot = plot, clean = clean, write = write)

}
