#' Title
#'
#' @param raw
#' @param ppm
#' @param snthresh
#' @param noise
#' @param prefilter
#' @param peakwidth
#' @param integrate
#' @param verbose
#' @param fitGauss
#' @param clean
#'
#' @return
#' @export
#'
#' @examples
pickPEAKS <- function(raw, ppm, snthresh, noise, prefilter, peakwidth, integrate, verbose, fitGauss, clean = TRUE) {

  stime <- Sys.time()

  ## set centwave parameters and find peaks
  CWParam <- xcms::CentWaveParam(ppm = ppm,
                                 snthresh = snthresh,
                                 noise = noise,
                                 prefilter = prefilter,
                                 peakwidth = peakwidth,
                                 integrate = integrate,
                                 verboseColumns = verbose,
                                 fitgauss = fitGauss)
  res <- xcms::findChromPeaks(object = raw, param = CWParam)
  pks <- data.frame(xcms::chromPeaks(res))

  ## order by intensity
  pks <- pks[order(pks$into, decreasing = T), ]

  ## give each peak an ID based on its intensity order
  pks$pid <- 1:nrow(pks)

  if(clean == TRUE) {

    ## clean-up unneccesary duplicating peaks
    message("'clean' set to TRUE. Removing duplicating peaks.")

    pks <- pks %>%
      group_by(rt, mz) %>%
      arrange(pid) %>%
      ## take only first of the two identical peaks which are in the same component
      filter(row_number()== 1) %>%
      ungroup() %>%
      ## update IDs
      mutate(pid = row_number()) %>%
      data.frame()

  }

  etime <- Sys.time()
  print(difftime(etime, stime))
  return(pks)
}


