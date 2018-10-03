context("groupPEAKS_paral")

test_that("Functions used by groupPEAKS_paral() returs correct output", {

  pickPEAKS_out <- pickPEAKS(raw = test_raw, fname = test_fname, cwt = cwt)
  extractEIC_out <- extractEIC(raw = test_raw, pks = test_pks_rd)

  ## (1) pickPEAKS()
  ## is the same number of peaks returned?
  expect_equal(nrow(pickPEAKS_out), nrow(test_pks_rd))
  ## are peaks in the table the same?
  expect_true(all(pickPEAKS_out$mz == test_pks_rd$mz))

  ## (2) extractEIC
  ## is the same number of EICs returned?
  expect_equal(length(extractEIC_out), length(test_eic_rd))

})





