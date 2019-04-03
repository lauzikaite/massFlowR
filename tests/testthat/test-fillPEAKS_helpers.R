context("helper functions for fillPEAKS method")

test_that("all helper functions for fillPEAKS work ", {
  # modelPEAKS --------------------------------------------------------------------------------------------------------
  ## use large study
  tmp <- buildTMP(file = large_meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  tmp <- validPEAKS(tmp, out_dir = data_dir, ncores = 2)
  
  ## model rt/mz values for a single peakgr that has missing values
  pcs_1 <- tmp@valid[tmp@valid$pcs == 1, "peakid"]
  data_2 <- tmp@data[[2]]
  peak_miss <- pcs_1[which(!pcs_1%in% data_2[data_2$tmp_peakgr ==1 ,"tmp_peakid" ])]
  peaks_ori <- tmp@peaks[[match(peak_miss, names(tmp@peaks))]]

  modelPEAKS_out <- modelPEAKS(p = peak_miss,
    vars = c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax"),
    object = tmp
    )
  expect_true(class(modelPEAKS_out) == "data.frame")
  expect_true(nrow(modelPEAKS_out) == nrow(tmp@samples))
  expect_equal(
    peaks_ori[!is.na(peaks_ori$mz), names(modelPEAKS_out)],
    modelPEAKS_out[!is.na(peaks_ori$mz), ]
    )
  expect_true(all(apply(modelPEAKS_out[is.na(peaks_ori$mz), ], 1:2, is.numeric)))
  
  # smoothVALUE --------------------------------------------------------------------------------------------------------
  
})



