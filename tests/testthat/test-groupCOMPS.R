context("groupCOMPS")

test_that("groupCOMPS of two identical tables returns correct template back", {

  getCOMPS(files = faahko_file, out_dir = massFlowR_dir, cwt = paramCWT)
  f <- list.files(massFlowR_dir, pattern = "_pks-comps.txt", full.names = T)
  tmp_with2 <- groupCOMPS(files = c(f, f), mz_err = 0.01, rt_err = 10)
  tmp_with3 <- groupCOMPS(files = c(f, f, f), mz_err = 0.01, rt_err = 10)

  ## expect that all peaks in one table will be grouped with themselves in the second table
  expect_true(all(tmp_with2$new_tmp$cid == tmp_with2$new_tmp$comp))
  expect_true(tmp_with2$new_tmp %>% filter(is.na(comp)) %>% nrow == 0)

  ## expect that grouping 3 files will still return identical tables
  expect_true(all(tmp_with3$new_tmp$cid == tmp_with3$new_tmp$comp))
  expect_true(tmp_with3$new_tmp %>% filter(is.na(comp)) %>% nrow == 0)

})
