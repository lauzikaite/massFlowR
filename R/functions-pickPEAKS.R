#' @title Peak detection using the centWave method
#' @description High resolution LC/MS datafiles are pick-peaked using centWave algorithm, implemented by \emph{xcms} package [Tautenhahn 2008].
#'
#' Required input arguments are identical to the ones required by \emph{CentWaveParam()} function in \emph{xcms} package.
#'
#' Object of class \code{XCMSnExp} is converted to data frame and un-nneccasary duplicating peaks are removed.
#'
#' @param raw \code{OnDiskMSnExp} class object.
#' @param ppm numeric defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.
#' @param snthresh numeric defining the signal to noise ratio cutoff.
#' @param noise numeric allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection).
#' @param prefilter numeric: c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I.
#' @param peakwidth numeric with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds.
#' @param integrate numeric for integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact.
#' @param fitGauss logical whether or not a Gaussian should be fitted to each peak.
#'
#' @return Function returns a data.frame class object with picked-peaks and centWave results.
#' @export
#'
#' @examples
#' raw <- MSnbase::readMSData(dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)[[1]], mode = "onDisk")
#' pks <- pickPEAKS(raw = raw, ppm = 25, snthresh = 10, noise = 0, prefilter = c(3, 100), peakwidth = c(30, 80), integrate = F, fitGauss = F)
#'

pickPEAKS <- function(raw,
                      ppm = xcms::CentWaveParam()@ppm,
                      snthresh = xcms::CentWaveParam()@snthresh,
                      noise = xcms::CentWaveParam()@noise,
                      prefilter = xcms::CentWaveParam()@prefilter,
                      peakwidth = xcms::CentWaveParam()@peakwidth,
                      integrate = xcms::CentWaveParam()@integrate,
                      fitGauss = xcms::CentWaveParam()@fitGauss) {

  ## set centwave parameters and find peaks
  CWParam <- xcms::CentWaveParam(ppm = ppm,
                                 snthresh = snthresh,
                                 noise = noise,
                                 prefilter = prefilter,
                                 peakwidth = peakwidth,
                                 integrate = integrate,
                                 verboseColumns = TRUE,
                                 fitgauss = fitGauss)
  res <- xcms::findChromPeaks(object = raw, param = CWParam)
  pks <- data.frame(xcms::chromPeaks(res))

  ## order by intensity
  pks <- pks[order(pks$into, decreasing = T), ]

  ## give each peak an ID based on its intensity order
  pks$pid <- 1:nrow(pks)

  pks_dup <- pks %>%
    group_by(rt, mz) %>%
    summarise(n = n()) %>%
    filter(n > 1) %>%
    nrow()

  message(pks_dup, "  duplicating peaks were removed.")

  pks <- pks %>%
    group_by(rt, mz) %>%
    arrange(pid) %>%
    ## take only first of the two identical peaks which are in the same component
    filter(row_number()== 1) %>%
    ungroup() %>%
    ## update IDs
    mutate(pid = row_number()) %>%
    data.frame()

  return(pks)
}


