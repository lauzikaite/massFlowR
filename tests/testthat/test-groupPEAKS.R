context("groupPEAKS")

test_that("Functions used by groupPEAKS_paral() returs correct output", {
  ## (1) readDATA
  ## is raw data read correctly
  readDATA_out <- readDATA(f = test_fname)
  expect_true(class(readDATA_out) == "OnDiskMSnExp")
  
  ## (2) pickPEAKS
  ## are the same peaks returned in the table?
  pickPEAKS_out <-
    pickPEAKS(raw = test_raw, fname = test_fname, cwt = cwt)
  expect_equal(nrow(pickPEAKS_out), nrow(test_pks_rd))
  expect_true(all(pickPEAKS_out$mz == test_pks_rd$mz))
  expect_true(all(pickPEAKS_out$peakid == test_pks_rd$peakid))
  
  ## (3) extractEIC
  ## is the same number of EICs returned?
  extractEIC_out <- extractEIC(raw = test_raw, pks = test_pks_rd)
  expect_equal(length(extractEIC_out), length(test_eic_rd))
  expect_equal(class(extractEIC_out[[1]])[1], "Chromatogram")
  
  ## (4) do_groupPEAKS
  do_groupPEAKS_out <-
    do_groupPEAKS(
      pks = test_pks_rd,
      eic = test_eic_rd,
      out_dir = data_dir,
      fname = test_basename,
      thr = 0.95,
      return = T
    )
  ## are only peakgrs with 2 or more peaks returned?
  expect_true(min(table(do_groupPEAKS_out$peakgr)) >= 2)
  ## are all peaks assigned to a peakgr
  expect_false(any(is.na(do_groupPEAKS_out$peakgr)))
  ## peakids are newly generated
  expect_equal(do_groupPEAKS_out$peakid, 1:nrow(do_groupPEAKS_out))
  
})

## functions specific to do_groupPEAKS
test_that("corEIC returns correct table", {
  ## get co-eluting peak indeces for peak no 228 (random)
  test_ind <-
    c(test_pks_rd[228, "scpos"] - 1, test_pks_rd[228, "scpos"] + 1)
  test_co <- test_pks_rd %>%
    filter(between(.data$scpos, test_ind[1], test_ind[2])) %>%
    pull(.data$peakid)
  test_cormat <- massFlowR:::getCORmat(ind = test_co)
  test_cormat_pair <- as.matrix(test_cormat)[1,]
  corEIC_out <-
    massFlowR:::corEIC(pair = test_cormat_pair, eic = test_eic_rd)
  expect_equal(round(corEIC_out, 7), 0.9061045)
  test_cormat$weight <- apply(test_cormat,
                              1,
                              FUN = massFlowR:::corEIC,
                              eic = test_eic_rd)
  expect_true(all(is.numeric(test_cormat$weight)))
  expect_true(all(test_cormat$weight >= 0))
})
