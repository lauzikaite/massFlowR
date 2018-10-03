context("groupPEAKSspec")

## get co-eluting peak indeces for peak no 228
test_pks_co_scans <- c(test_pks_rd[228,"scpos"] - 1, test_pks_rd[228, "scpos"] + 1)
test_pks_co <- test_pks_rd %>%
  filter(between(.data$scpos, test_pks_co_scans[1], test_pks_co_scans[2])) %>%
  pull(.data$peakid)

## (1) test that builCOR() correctly calculates correlation between co-eluting peaks
test_that("buildCOR returns correct table", {

  buildCOR_out <- buildCOR(co_ind = test_pks_co, eic = test_eic_rd, pearson = TRUE)

  ## all co-eluting peaks are in the output?
  expect_true(all(buildCOR_out$x %in% test_pks_co, buildCOR_out$y %in% test_pks_co))
  expect_true(is.numeric(buildCOR_out$cor))

})
