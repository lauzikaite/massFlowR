#' @title Obtain extracted ion chromatograms for all picked peaks
#' @description Obtain extracted ion chromatograms (EIC) for all picked peaks, provided by \emph{pickPEAKS} function.
#'
#' EICs are stored in a list of length equal to the number of peaks. Individual peaks' EICs can be plotted using base graphics.
#'
#' @param raw \code{OnDiskMSnExp} class object for the single datafile of interest.
#' @param pks \code{DataFrame} object generated with picked peaks.
#'
#' @return Function returns a \code{list} with an EIC for eack peak in the \code{pks} table.
#'
#' @examples
#' raw <- MSnbase::readMSData(dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)[[1]], mode = "onDisk")
#' pks <- pickPEAKS(raw = raw, ppm = 25, snthresh = 10, noise = 0, prefilter = c(3, 100), peakwidth = c(30, 80), integrate = F, fitGauss = F)
#' eic <- extractEIC(raw = raw, pks = pks)
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
