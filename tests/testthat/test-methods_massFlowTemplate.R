context("methods for massFlowTemplate")

# checkNEXT ---------------------------------------------------------------------------------------------------
test_that("checkNEXT() returns the correct filename", {
  
  tmp <- buildTMP(file = experiment_file)
  expected_next <- grouped_files[2]
  expect_equal(expected_next, checkNEXT(tmp))
  
})


# align identical tables---------------------------------------------------------------------------------------------------
test_that("alignment of two identical samples via alignPEAKS() is correct", {
  
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

test_that("alignment of two different samples via alignPEAKS() is correct", {
  
  tmp <- buildTMP(file = experiment_file, rt_err = 10)
  tmp <- alignPEAKS(tmp)
  tmp_tmp <- tmp@tmp
  tmp_data1 <- tmp@data[[1]]
  tmp_data2 <- tmp@data[[2]]
  
  # all samples were aligned
  expect_true(all(tmp@samples[,"aligned"])) 
  
  ####---- returned template is correct
  # all template peaks come from the aligned samples
  expect_true(all(tmp_data1$tmp_peakid %in% tmp_tmp$peakid))
  expect_true(all(tmp_data2$tmp_peakid %in% tmp_tmp$peakid))
  
  # for peaks detected ONLY in last sample:
  # intensity values in the final template equals intensity values in the latest sample
  peakid_last <- match(tmp_data2$tmp_peakid, tmp_tmp$peakid)
  expect_true(all(tmp_tmp[peakid_last,"into"] == tmp_data2$into))
  
  # for peaks detected ONLY in first sample:
  # intensity values in the final template equals intensity values in the first sample 
  peakid_first <- tmp_tmp$peakid[match(tmp_data1$tmp_peakid, tmp_tmp$peakid)] 
  peakid_first <- peakid_first[which(!peakid_first %in% tmp_data2$tmp_peakid)] 
  peakid_first_id <- match(peakid_first, tmp_tmp$peakid)
  expect_true(all(tmp_tmp[peakid_first_id,"into"] == tmp_data1$into[which(tmp_data1$tmp_peakid %in%  peakid_first)]))
  
  ####---- returned doi files are correct
  ## all peaks are returned
  data1 <- read.csv(grouped_files[1], header = T, stringsAsFactors = F)
  data2 <- read.csv(grouped_files[2], header = T, stringsAsFactors = F)
  
  expect_true(all(data1$peakid %in% tmp_data1$peakid))
  expect_equal(nrow(data1), nrow(tmp_data1))
  expect_true(all(data2$peakid %in% tmp_data2$peakid))
  expect_equal(nrow(data2), nrow(tmp_data2))
  
  ## outpout must have slightly different rt values
  ## is the new_rt the average between tmp and data2?
  for (p in tmp_data2$tmp_peakid) {
    if (p %in% tmp_data1$tmp_peakid) {
      rt_1 <- tmp_data1[which(tmp_data1$tmp_peakid == p),"rt"]
      rt_1_new <- tmp_data1[which(tmp_data1$tmp_peakid == p),"new_rt"]
      expect_equal(rt_1, rt_1_new)
      rt_2 <- tmp_data2[which(tmp_data2$tmp_peakid == p),"rt"]
      rt_2_new <- tmp_data2[which(tmp_data2$tmp_peakid == p),"new_rt"]
      rt_tmp <- tmp_tmp[which(tmp_tmp$peakid == p),"rt"]
      expect_equal(rt_2_new, rt_tmp)
      expect_equal(rt_2_new, median(rt_1, rt_2))
    }
  }
  
})
  
# align almost identical tables---------------------------------------------------------------------------------------------------
test_that("alignment of two ALMOST identical samples via alignPEAKS() is correct", {
  
  tmp <- buildTMP(file = experiment_mess_file, rt_err = 10)
  tmp <- alignPEAKS(tmp)
  
  tmp_tmp <- tmp@tmp
  tmp_data1 <- tmp@data[[1]]
  tmp_data2 <- tmp@data[[2]]
  
  ## peakgroup equal to "biggest_pkg" was present both in data1 and data2
  ## data2 has two copies of this peakgrop:
  # 1. is the original but with fewer peaks (n/2 - 1)
  # 2. is modifed with altered rt/into and fewer peaks (n/2)
  
  ## double checking whether the peakgroup-of-interest is what I expect
  pkg_grouped_1 <- tmp_data1[which(tmp_data1$tmp_peakgr == biggest_pkg),]
  pkg_grouped_2 <- tmp_data2[which(tmp_data2$tmp_peakgr == biggest_pkg),]
  expect_equal(nrow(pkg_grouped_2), (nrow(pkg_grouped_1)/2 - 1))
  expect_true(all(pkg_grouped_2$tmp_peakid %in% pkg_grouped_1$tmp_peakid))
  
  ## checking whether addDOI correctly "kicked-out" previously grouped messy peakgr
  expect_true(all(is.na(tmp_data2[which(tmp_data2$tmp_peakgr == messy_pkg),"cos"])))
  expect_false(any(tmp_tmp[which(tmp_tmp$peakgr == messy_pkg),"peakid"] %in%  tmp_data1$tmp_peakid))
})

