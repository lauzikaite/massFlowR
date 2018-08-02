context("groupCOMPS")

test_that("groupCOMPS of two identical tables returns correct template back", {
  paramCWT <- xcms::CentWaveParam(ppm = ppm,
                                  snthresh = snthresh,
                                  noise = noise,
                                  prefilter = prefilter,
                                  peakwidth = peakwidth,
                                  integrate = integrate,
                                  fitgauss = fitGauss,
                                  verboseColumns = TRUE)
  getCOMPS(files = fl[1], out_dir = out_dir, cwt = paramCWT)
  f <- list.files(out_dir, pattern = "_pks-comps-cls.txt", full.names = T)
  tmp <- groupCOMPS(files = c(f, f), mz_err = mz_err, rt_err = rt_err)

  ## expect that all peaks in one table will be grouped with themselves in the second table
  expect_true(all(tmp$new_tmp$cid == tmp$new_tmp$comp))
  expect_true(tmp$new_tmp %>% filter(is.na(comp)) %>% nrow == 0)

})
