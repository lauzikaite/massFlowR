#' @title Peak detection using the centWave method
#' @description High resolution, centroid LC/MS datafiles are pick-peaked using centWave algorithm, implemented by \emph{xcms} package [Tautenhahn 2008].
#'
#' Required input arguments are identical to the ones required by \emph{CentWaveParam} function in \emph{xcms} package.
#'
#' Object of class \code{XCMSnExp} is converted to \code{DataFrame} and duplicating peaks with the same \code{rt} and \code{mz} are merged into one.
#'
#' \emph{pickPEAKS} is applied for a single datafile.
#'
#' @param cwt \code{CentWaveParam} class object with parameters for peak-picking.
#' @param raw \code{OnDiskMSnExp} class object for the single datafile of interest.
#' @param fname \code{character} object specifying datafile name.
#' @param out_dir \code{character} object specifying directory where output data will be saved.
#'
#' @return Function returns a \code{DataFrame} object with picked-peaks and centWave details.
#' @export
#'
#' @examples

pickPEAKS <- function(raw, cwt, fname, out_dir, write = TRUE) {

  ## set centwave parameters and find peaks
  if (missing(cwt)) { stop("'cwt' has to be specified!") }
  if (class(cwt) != "CentWaveParam") { stop("'cwt' has to be 'CentWaveParam' object!") }

  if(write == TRUE) { if(missing(fname)) { stop("'fname' must be specified!") } }
  if(write == TRUE) { if(missing(out_dir)) { stop("'out_dir' must be specified!") } }

  res <- xcms::findChromPeaks(object = raw, param = cwt)
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

  if (write == T) {
    message(fname, " . Writing peaks table to txt file ...")
    write.table(pks, file = paste0(out_dir, fname, "_pks.txt"), col.names = T, quote = F, sep = "\t", row.names = F)
  }

  return(pks)

}


