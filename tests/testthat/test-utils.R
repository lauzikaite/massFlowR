context("utils functions")

test_that("cleanPEAKS removes duplicated peaks correctly", {
  
  ## removes duplicates in mz/rt (as in pickPEAKS)
  # use un-cleaned, original centwave output
  pks_unique <- unique(test_pks[,c("mz", "rt")])
  pks_clean <- lapply(1:nrow(pks_unique),
                      FUN = cleanPEAKS,
                      dt_unique = pks_unique,
                      dt = test_pks)
  pks_clean <- do.call("rbindCLEAN", pks_clean)
  expect_true(all_equal(test_pks_rd %>% select(-peakid),
                        pks_clean))
  
  ## removes duplicated peakids (as in addDOI)
  # create peak table with duplicated peakids
  test_pks_dup <- test_pks_rd %>%
    bind_rows(., test_pks_rd %>% filter(peakid %in% c(1,10,20)))
  # cleanup peak table
  pks_unique <- unique(test_pks_dup[,c("mz", "rt")])
  pks_clean <- lapply(1:nrow(pks_unique),
                      FUN = cleanPEAKS,
                      dt_unique = pks_unique,
                      dt = test_pks_dup)
  pks_clean <- do.call("rbindCLEAN", pks_clean)
  expect_true(any(duplicated(test_pks_dup$peakid)))
  expect_true(all_equal(test_pks_rd, pks_clean))
})

test_that("rbindCLEAN", {
    
  ## make a list of tables
  # replicate same table ntimes
  ntimes <- 10
  list_of_pks <- lapply(1:ntimes, function(dt) {
    test_pks_rd
  })
  dt_of_pks <- do.call("rbindCLEAN", list_of_pks)
  expect_true(all(rownames(dt_of_pks) == 1:nrow(dt_of_pks)))
  expect_equal(length(unique(dt_of_pks$peakid)), max(test_pks_rd$peakid))
  expect_equal(length(which(dt_of_pks$peakid == 1)), ntimes)
})



