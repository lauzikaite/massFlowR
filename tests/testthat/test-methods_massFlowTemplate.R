context("methods for massFlowTemplate")

test_that("checkNEXT returns the correct filename", {
  
  tmp <- buildTMP(file = experiment_file)
  expected_next <- grouped_files[2]
  expect_equal(expected_next, checkNEXT(tmp))
  
})


test_that("align two identical samples, without DB", {
  
  tmp <- buildTMP(file = experiment_dup_file, rt_err = 10)
  tmp <- alignPEAKS(tmp)
  tmp_data1 <- tmp@data[[1]]
  tmp_data2 <- tmp@data[[2]]
  
  # all samples were aligned
  expect_true(all(tmp@samples[,"aligned"])) 
  # aligned tables for identical samples are also identical
  expect_true(all_equal(
    tmp_data1 %>% select(-cos),
    tmp_data2 %>% select(-cos),
    # table 1 will have peakids as integers, table 2 will have peakids as numerics because addDOI converts to numeric
    convert = T))
  # generated tmp is identical to original table
  expect_equal(nrow(tmp@tmp), nrow(single_table))
  expect_true(all_equal(tmp@tmp[which(tmp@tmp$peakid %in% single_table$peakid),c("peakid", "mz", "rt", "into", "peakgr")],
                        single_table[,c("peakid", "mz", "rt", "into", "peakgr")]))

})

test_that("align two different samples, without DB", {
  
  tmp <- buildTMP(file = experiment_file, rt_err = 10)
  tmp <- alignPEAKS(tmp)
  tmp_data1 <- tmp@data[[1]]
  tmp_data2 <- tmp@data[[2]]
  
  # all samples were aligned
  expect_true(all(tmp@samples[,"aligned"])) 
  
  # all template peaks come from the aligned samples
  expect_true(all(tmp_data1$tmp_peakid %in% tmp@tmp$peakid))
  expect_true(all(tmp_data2$tmp_peakid %in% tmp@tmp$peakid))
  
  # for peaks detected ONLY in last sample:
  # intensity values in the final template equals intensity values in the latest sample
  peakid_last <- match(tmp_data2$tmp_peakid, tmp@tmp$peakid)
  expect_true(all(tmp@tmp[peakid_last,"into"] == tmp_data2$into))
  
  # for peaks detected ONLY in first sample:
  # intensity values in the final template equals intensity values in the first sample 
  peakid_first <- tmp@tmp$peakid[match(tmp_data1$tmp_peakid, tmp@tmp$peakid)] 
  peakid_first <- peakid_first[which(!peakid_first %in% tmp_data2$tmp_peakid)] 
  peakid_first_id <- match(peakid_first, tmp@tmp$peakid)
  expect_true(all(tmp@tmp[peakid_first_id,"into"] == tmp_data1$into[which(tmp_data1$tmp_peakid %in%  peakid_first)]))
  
})

  