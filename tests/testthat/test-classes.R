context("data classes and classes constructors")

# buildTMP ------------------------------------------------------------------------------------------------------
test_that("initialise massFlowTemplate class object using default parameters", {
  tmp <- buildTMP(file = meta_fname, out_dir = data_dir)

  ## Check object slots
  expect_true(class(tmp) == "massFlowTemplate")

  ## Check object slot values
  expect_true(all(tmp@samples$raw_filepath == metadata$raw_filepath))
  expect_true(tmp@samples[which(tmp@samples$aligned == TRUE), "proc_filepath"] == metadata$proc_filepath[1])
  expect_true(all(tmp@tmp[match(tmp@tmp$peakid, single_table$peakid), "mz"] == single_table$mz))
  expect_true(all(tmp@tmp[match(tmp@tmp$peakid, single_table$peakid), "into"] == single_table$into))

  data1 <- tmp@data[[1]]
  expect_true(all(data1$peakid %in% single_table$peakid))
  expect_true(all(data1$peakid[match(data1$peakid, single_table$peakid)] == single_table$peakid))
  expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == 2, tmp@params$bins == 0.05))
})

# loadALIGNED------------------------------------------------------------------------------------------------------
test_that("create massFlowTemplate class object and fill it with aligned samples", {
  tmp <- buildTMP(file = meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  tmp <- loadALIGNED(file = file.path(data_dir, "aligned.csv"),
                     template = file.path(data_dir, "template.csv"),
                     rt_err = rt_err)
  
  expect_false(peaksVALIDATED(tmp))
  expect_true(validObject(tmp))
  expect_equal(rt_err, tmp@params$rt_err)
}) 

