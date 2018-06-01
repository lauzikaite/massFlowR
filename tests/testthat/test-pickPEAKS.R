context("pickPEAKS")

paramCWT <- xcms::CentWaveParam(ppm = ppm,
                                snthresh = snthresh,
                                noise = noise,
                                prefilter = prefilter,
                                peakwidth = peakwidth,
                                integrate = integrate,
                                fitgauss = fitGauss,
                                verboseColumns = TRUE)

## get full peak table
get_full_pks <- function(files, cwt) {

  raw <-  MSnbase::readMSData(files = files, mode = "onDisk")
  res <- xcms::findChromPeaks(object = raw, param = cwt)
  pks <- data.frame(xcms::chromPeaks(res))
  return(pks)
}

## get peak table without dupolicating peaks
get_reduced_pks <- function(files, cwt){

  ## get peak table using basic xcms functionality
  pks <- get_full_pks(files = files,
                      cwt = cwt)
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



test_that("pickPEAKS returns correct table", {

  ## select one file
  f <- files[1]

  ## test with defined helper functions
  out <- get_reduced_pks(files = f, cwt = paramCWT)

  ## test with package pickPEAKS() function
  raw <-  MSnbase::readMSData(files = f, mode = "onDisk")
  pickPEAKS_out <- pickPEAKS(raw = raw, cwt = paramCWT, write = FALSE)

  expect_equal(nrow(pickPEAKS_out), nrow(out))

})
