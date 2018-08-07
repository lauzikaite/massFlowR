context("groupCOMPS")

test_that("groupCOMPS of two identical tables returns correct template back", {

  getCOMPS(files = faahko_file, out_dir = massflowR_dir, cwt = paramCWT)
  f <- list.files(massflowR_dir, pattern = "_pks-comps.txt", full.names = T)
  tmp <- groupCOMPS(files = c(f, f), mz_err = 0.01, rt_err = 10)

  ## expect that all peaks in one table will be grouped with themselves in the second table
  expect_true(all(tmp$new_tmp$cid == tmp$new_tmp$comp))
  expect_true(tmp$new_tmp %>% filter(is.na(comp)) %>% nrow == 0)

})
