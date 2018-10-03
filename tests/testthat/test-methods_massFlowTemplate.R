context("methods for massFlowTemplate")

test_that("align two identical samples, without DB", {
  
  tmp <- buildTMP(file = experiment_dup_file, rt_err = 10)
  tmp <- alignSAMPLES(tmp)
  tmp_data1 <- tmp@data[[1]]
  tmp_data2 <- tmp@data[[2]]
  
  # all samples were aligned
  expect_true(all(tmp@samples[,"aligned"])) 
  # aligned tables for identical samples are also identical
  expect_true(all_equal(tmp_data1 %>% select(-cos),
                        tmp_data2 %>% select(-cos)))
  # generated tmp is identical to original table
  expect_equal(nrow(tmp@tmp), nrow(single_table))
  expect_true(all_equal(tmp@tmp[which(tmp@tmp$peakid %in% single_table$peakid),c("peakid", "mz", "rt", "into", "peakgr")],
                        single_table[,c("peakid", "mz", "rt", "into", "peakgr")]))

})
  