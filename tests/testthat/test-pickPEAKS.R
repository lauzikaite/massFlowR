context("pickPEAKS")

test_that("pickPEAKS returns correct table", {
  pks_reduced <- get_reduced_pks(files = files,
                                 ppm = ppm,
                                 snthresh = snthresh,
                                 noise = noise,
                                 prefilter = prefilter,
                                 peakwidth = peakwidth,
                                 integrate = integrate,
                                 verbose = verbose,
                                 fitGauss = fitGauss)
  raw <-  MSnbase::readMSData(files = files, mode = "onDisk")
  pickPEAKS_out <- pickPEAKS(raw = raw,
                             ppm = ppm,
                             snthresh = snthresh,
                             noise = noise,
                             prefilter = prefilter,
                             peakwidth = peakwidth,
                             integrate = integrate,
                             verbose = verbose,
                             fitGauss = fitGauss,
                             clean == clean)

  expect_equal(nrow(pickPEAKS_out), nrow(pks_reduced))

})
