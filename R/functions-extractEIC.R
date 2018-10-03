#' @title Obtain extracted ion chromatograms for all picked peaks
#' 
#' @description Obtain extracted ion chromatograms (EIC) for all picked peaks, provided by \emph{pickPEAKS} function.
#' 
#' @param raw \code{OnDiskMSnExp} class object for the single datafile of interest.
#' @param pks \code{DataFrame} object generated with picked peaks.
#'
#' @return Function returns a \code{list} with an EIC for eack peak in the \code{pks} table.
#' Individual peaks' EICs can be plotted using base graphics.
#'
#' @examples
#' fname <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE)[[1]]
#' raw <- MSnbase::readMSData(fname, mode = "onDisk")
#' cwt <-  xcms::CentWaveParam(ppm = 25, snthresh = 10, noise = 0,
#' prefilter = c(3, 100), peakwidth = c(30, 80))
#' pks <- pickPEAKS(raw = raw, fname = basename(fname), cwt = cwt)
#' eic <- extractEIC(raw = raw, pks = pks)
#' 
#' ## now you can plot a single peak's EIC
#' plot(eic[[1]])

extractEIC <- function(raw, pks) {

  eic <- xcms::chromatogram(raw,
                            rt = data.frame(
                              rt_lower = pks$rtmin,
                              rt_upper = pks$rtmax),
                            mz = data.frame(
                              mz_lower = pks$mzmin,
                              mz_upper = pks$mzmax))
  clean_eic <- lapply(1:nrow(eic), function(ch) {
    clean(eic[ch, ], na.rm = T)
  })
  return(clean_eic)
}
