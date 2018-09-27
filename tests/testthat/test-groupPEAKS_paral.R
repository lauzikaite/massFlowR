context("groupPEAKS_paral")

test_that("Functions used by groupPEAKS_paral() returs correct output", {

  pickPEAKS_out <- pickPEAKS(raw = faahko_raw, fname = faahko_fname, cwt = cwt, out_dir = massFlowR_dir)
  extractEIC_out <- extractEIC(raw = faahko_raw, pks = faahko_pks_rd)

  ## (1) pickPEAKS()
  ## is the same number of peaks returned?
  expect_equal(nrow(pickPEAKS_out), nrow(faahko_pks_rd))
  ## are peaks in the table the same?
  expect_true(all(pickPEAKS_out$mz == faahko_pks_rd$mz))

  ## (2) extractEIC
  ## is the same number of EICs returned?
  expect_equal(length(extractEIC_out), length(faahko_eic_rd))



})





