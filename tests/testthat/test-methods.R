context("methods for massFlowTemplate")

# basic methods -----------------------------------------------------------------------------------------------------
test_that("basic methods for massFlowTemplate", {
  tmp <- buildTMP(file = meta_fname, out_dir = data_dir)
  expect_true(validmassFlowTemplate(tmp))
  expect_true(filepath(tmp) == meta_fname)
  expect_false(peaksVALIDATED(tmp))
})


# checkNEXT ---------------------------------------------------------------------------------------------------
test_that("checkNEXT() returns the correct filename", {
  tmp <- buildTMP(file = meta_fname, out_dir = data_dir)
  expected_next <- grouped_fnames[2]
  expect_equal(expected_next, checkNEXT(tmp))
})

# alignPEAKS identical tables---------------------------------------------------------------------------------------------------
test_that("alignment of two identical samples via alignPEAKS() is correct", {
  tmp <- buildTMP(file = dup_meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  tmp_data1 <- tmp@data[[1]]
  tmp_data2 <- tmp@data[[2]]
  
  # all samples were aligned
  expect_true(all(tmp@samples[,"aligned"])) 
  # aligned tables for identical samples are also identical
  expect_true(all(tmp_data1[match(tmp_data1$peakid, tmp_data2$peakid), "mz"] == tmp_data2$mz))
  expect_true(all(tmp_data1[match(tmp_data1$peakid, tmp_data2$peakid), "rt"] == tmp_data2$rt))
  
  # generated tmp is identical to original table
  expect_equal(nrow(tmp@tmp), nrow(single_table))
  expect_true(all(tmp@tmp[match(tmp@tmp$peakid, single_table$peakid), c("peakid", "mz", "rt", "into", "peakgr")] ==
                    single_table[,c("peakid", "mz", "rt", "into", "peakgr")]))
})

# alignPEAKS different tables---------------------------------------------------------------------------------------------------
test_that("alignment of two different samples via alignPEAKS() is correct", {
  tmp <- buildTMP(file = meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
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
  data1 <- read.csv(grouped_fnames[1], header = T, stringsAsFactors = F)
  data2 <- read.csv(grouped_fnames[2], header = T, stringsAsFactors = F)
  
  expect_true(all(data1$peakid %in% tmp_data1$peakid))
  expect_equal(nrow(data1), nrow(tmp_data1))
  expect_true(all(data2$peakid %in% tmp_data2$peakid))
  expect_equal(nrow(data2), nrow(tmp_data2))
})
  
# alignPEAKS almost identical tables---------------------------------------------------------------------------------------------------
test_that("alignment of two ALMOST identical samples via alignPEAKS() is correct", {
  tmp <- buildTMP(file = noisy_meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  
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
  
  ## checking whether addDOI correctly "kicked-out" previously grouped noisy peakgr
  expect_true(all(is.na(tmp_data2[which(tmp_data2$tmp_peakgr == noisy_pkg), "cos"])))
  expect_false(any(tmp_tmp[which(tmp_tmp$peakgr == noisy_pkg), "peakid"] %in% tmp_data1$tmp_peakid))
})

# validPEAKS---------------------------------------------------------------------------------------------------
test_that("validPEAKS", {
  ####---- object with 2 samples will stop due to size
  tmp <- buildTMP(file = meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  expect_error(validPEAKS(tmp, out_dir = data_dir, ncores = 2),
               "object has 2 samples\n minimum 3 samples are required for validation.")
  
  ####---- use large study
  tmp <- buildTMP(file = large_meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  validPEAKS_out <- validPEAKS(tmp, out_dir = data_dir, ncores = 2)
  expect_true(class(validPEAKS_out) == "massFlowTemplate")
  expect_true(peaksVALIDATED(validPEAKS_out))
  # check object slots
  expect_true(nrow(validPEAKS_out@valid) > 0)
  expect_true(length(validPEAKS_out@peaks) == nrow(validPEAKS_out@valid)) # an entry for every peak in the final template
  expect_true(length(validPEAKS_out@values) == nrow(validPEAKS_out@samples)) # an entry for every sample
  expect_true(NA %in% unlist(lapply(validPEAKS_out@values, "[[", "into"))) # some missed peaks with NA
  expect_true(all(sapply(validPEAKS_out@peaks, nrow) == nrow(tmp@samples))) # a row for every peak for every sample
  # were files written?
  expect_true(length(grep("intensity_data.csv", list.files(data_dir))) > 0)
  expect_true(length(grep("peaks_data.csv", list.files(data_dir))) > 0)
  expect_true(length(grep("sample_data.csv", list.files(data_dir))) > 0)
  expect_true(length(grep("object.RDS", list.files(data_dir))) > 0)
})

# fillPEAKS---------------------------------------------------------------------------------------------------
test_that("fillPEAKS", {
  ####---- use large study
  tmp <- buildTMP(file = large_meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  expect_error(fillPEAKS(tmp)) # needs validation first
  tmp <- validPEAKS(tmp, out_dir = data_dir)
  fillPEAKS_out <- fillPEAKS(tmp, out_dir = data_dir)
  expect_true(class(fillPEAKS_out) == "massFlowTemplate")
  expect_true(peaksVALIDATED(fillPEAKS_out))
  expect_true(length(fillPEAKS_out@values) == nrow(fillPEAKS_out@samples))
  expect_true(is.numeric(unlist(lapply(fillPEAKS_out@values, "[[", "into"))))
  expect_false(NA %in% unlist(lapply(fillPEAKS_out@values, "[[", "into"))) # NO NAs
  expect_true("into"  %in% colnames(fillPEAKS_out@valid)) # new column added with median into values
  # were files written?
  expect_true(length(grep("filled_intensity_data.csv", list.files(data_dir))) > 0)
  expect_true(length(grep("final_peaks_data.csv", list.files(data_dir))) > 0)
  expect_true(length(grep("sample_data.csv", list.files(data_dir))) > 0)
  expect_true(length(grep("object.RDS", list.files(data_dir))) > 0)
  
})
