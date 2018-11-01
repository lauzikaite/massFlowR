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
#' @export
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
      message("pickPEAKS fail. failing mzML file: ", fname)
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



