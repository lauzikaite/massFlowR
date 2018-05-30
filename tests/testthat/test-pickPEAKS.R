context("pickPEAKS")

## get peak table without dupolicating peaks
get_reduced_pks <- function(files = files,
                            ppm = ppm,
                            snthresh = snthresh,
                            noise = noise,
                            prefilter = prefilter,
                            peakwidth = peakwidth,
                            integrate = integrate,
                            fitGauss = fitGauss){

  ## get peak table using basic xcms functionality
  pks <- get_full_pks(files = files,
                      ppm = ppm,
                      snthresh = snthresh,
                      noise = noise,
                      prefilter = prefilter,
                      peakwidth = peakwidth,
                      integrate = integrate,
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



test_that("pickPEAKS returns correct table", {

  out <- get_reduced_pks(files = files,
                                 ppm = ppm,
                                 snthresh = snthresh,
                                 noise = noise,
                                 prefilter = prefilter,
                                 peakwidth = peakwidth,
                                 integrate = integrate,
                                 fitGauss = fitGauss)
  ## run pickPEAKS() function
  raw <-  MSnbase::readMSData(files = files, mode = "onDisk")
  pickPEAKS_out <- pickPEAKS(raw = raw,
                             ppm = ppm,
                             snthresh = snthresh,
                             noise = noise,
                             prefilter = prefilter,
                             peakwidth = peakwidth,
                             integrate = integrate,
                             fitGauss = fitGauss)

  expect_equal(nrow(pickPEAKS_out), nrow(out))

})
