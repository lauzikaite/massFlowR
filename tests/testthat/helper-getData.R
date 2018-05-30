## build components using faahKO package data

files <- dir(system.file("cdf", package = "faahKO"),
            full.names = TRUE, recursive = TRUE)[[1]]
noise <- 1000
peakwidth <-  c(30, 80)
prefilter <- c(3, 100)
snthresh <- 10
ppm <- 25
integrate <- 1
fitGauss <- FALSE
match <- 1
Pearson <- TRUE

get_full_pks <- function(files, ppm, snthresh, noise, prefilter, peakwidth, integrate, fitGauss) {

  raw <-  MSnbase::readMSData(files = files, mode = "onDisk")
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
  return(pks)
}


