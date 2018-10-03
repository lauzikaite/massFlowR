#' @aliases pickPEAKS
#' 
#' @title Peak detection using the centWave method
#' 
#' @description High resolution, centroid LC/MS datafiles are pick-peaked using centWave algorithm, implemented by \emph{xcms} package [Tautenhahn 2008].
#' 
#' \emph{pickPEAKS} is applied to a single datafile.
#'
#' @param cwt \code{CentWaveParam} class object with parameters for peak-picking. Object can be created by the \emph{xcms::CentWaveParam} function.
#' @param raw \code{OnDiskMSnExp} class object for the single LC-MS spectrum of interest. Such object can be created by the \emph{MSnbase::readMSData} function.
#' @param fname \code{character} object specifying datafile name.
#'
#' @return Function returns a \code{DataFrame} object with picked-peaks and centWave details.
#' 
#' @export
#'

pickPEAKS <- function(raw, cwt, fname) {

  message("Peak-picking file: ", fname, " ...")

  res <- xcms::findChromPeaks(object = raw, param = cwt)
  pks <- data.frame(xcms::chromPeaks(res))

  ## order by intensity
  pks <- pks[order(pks$into, decreasing = T), ]

  ## assign each peak with id, based on its intensity order
  pks$peakid <- 1:nrow(pks)

  ## remove artefactural, duplicating peaks
  pks <- pks %>%
    group_by(.data$rt, .data$mz) %>%
    arrange(.data$peakid) %>%
    ## take only first of the two identical peaks
    filter(row_number() == 1) %>%
    ungroup() %>%
    ## update peak id
    mutate(peakid = row_number()) %>%
    data.frame()

  return(pks)

}


