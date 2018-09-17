#' @title Peak detection using the centWave method
#' @description High resolution, centroid LC/MS datafiles are pick-peaked using centWave algorithm, implemented by \emph{xcms} package [Tautenhahn 2008].
#'
#' Required input arguments are identical to the ones required by \emph{CentWaveParam} function in \emph{xcms} package.
#'
#' Object of class \code{XCMSnExp} is converted to \code{DataFrame} and duplicating peaks with the same \code{rt} and \code{mz} are merged into one.
#'
#' \emph{pickPEAKS} is applied for a single datafile.
#'
#' @param cwt \code{CentWaveParam} class object with parameters for peak-picking. Object can be created by the \emph{xcms::CentWaveParam} function.
#' @param raw \code{OnDiskMSnExp} class object for the single LC-MS spectrum of interest. Such object can be created by the \emph{MSnbase::readMSData} function.
#' @param fname \code{character} object specifying datafile name.
#' @param out_dir \code{character} object specifying directory where output data will be saved.
#'
#' @return Function returns a \code{DataFrame} object with picked-peaks and centWave details.
#' @export
#'

pickPEAKS <- function(raw, cwt, fname, out_dir) {

  message("Peak-picking file: ", fname, " ...")

  res <- xcms::findChromPeaks(object = raw, param = cwt)
  pks <- data.frame(xcms::chromPeaks(res))

  ## order by intensity
  pks <- pks[order(pks$into, decreasing = T), ]

  ## assign each peak with id, based on its intensity order
  pks$peakid <- 1:nrow(pks)

  ## remove artefactural, duplicating peaks
  pks <- pks %>%
    group_by(rt, mz) %>%
    arrange(peakid) %>%
    ## take only first of the two identical peaks
    filter(row_number()== 1) %>%
    ungroup() %>%
    ## update peak id
    mutate(peakid = row_number()) %>%
    data.frame()

  write.csv(pks, file = paste0(out_dir, "/", fname, "_peaks.csv"), quote = F, row.names = F)

  return(pks)

}


