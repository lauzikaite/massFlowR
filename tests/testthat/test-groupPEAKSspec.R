context("groupPEAKSspec")

## get co-eluting peak indeces for peak no 228
faahko_pks_co_scans <- c(faahko_pks_rd[228,"scpos"] - 1, faahko_pks_rd[228, "scpos"] + 1)
faahko_pks_co <- faahko_pks_rd %>%
  filter(between(scpos, faahko_pks_co_scans[1], faahko_pks_co_scans[2])) %>%
  pull(peakid)

## test that internal builCOR function correctly calculates correlation between co-eluting peaks

test_that("buildCOR returns correct table", {

  buildCOR_out <- buildCOR(co_ind = faahko_pks_co, eic = faahko_eic_rd, pearson = TRUE)

  ## all co-eluting peaks are in the output?
  expect_true(all(buildCOR_out$x %in% faahko_pks_co, buildCOR_out$y %in% faahko_pks_co))
  expect_true(is.numeric(buildCOR_out$cor))

})
