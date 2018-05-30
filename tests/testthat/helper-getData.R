## helper function to build components using faahKO package data

files <- dir(system.file("cdf", package = "faahKO"),
            full.names = TRUE, recursive = TRUE)[[1]]
noise <- 1000
peakwidth <-  c(30, 80)
prefilter <- xcms::CentWaveParam()@prefilter
snthresh <- xcms::CentWaveParam()@snthresh
ppm <- xcms::CentWaveParam()@ppm
integrate <- 1
verbose <- TRUE
fitGauss <- FALSE
match <- 1
Pearson <- TRUE
clean <- TRUE

get_full_pks <- function(files, ppm, snthresh, noise, prefilter, peakwidth, integrate, verbose, fitGauss) {

  raw <-  MSnbase::readMSData(files = files, mode = "onDisk")
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
  return(pks)
}

get_reduced_pks <- function(files = files,
                            ppm = ppm,
                            snthresh = snthresh,
                            noise = noise,
                            prefilter = prefilter,
                            peakwidth = peakwidth,
                            integrate = integrate,
                            verbose = verbose,
                            fitGauss = fitGauss){

  pks <- get_full_pks(files = files,
                      ppm = ppm,
                      snthresh = snthresh,
                      noise = noise,
                      prefilter = prefilter,
                      peakwidth = peakwidth,
                      integrate = integrate,
                      verbose = verbose,
                      fitGauss = fitGauss)

  pks <- pks[order(pks$into, decreasing = T), ]
  pks$pid <- 1:nrow(pks)

  pks_reduced <- pks %>%
    group_by(rt, mz) %>%
    arrange(pid) %>%
    filter(row_number()== 1) %>%
    ungroup() %>%
    mutate(pid = row_number()) %>%
    data.frame()

  return(pks_reduced)
}
