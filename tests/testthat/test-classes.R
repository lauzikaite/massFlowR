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
  ## check all params values
  expect_true(length(tmp@params) == 5)
  expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == 2, tmp@params$bins == 0.05, tmp@params$cutoff == 0, tmp@params$realtime == FALSE))
})

# loadALIGNED------------------------------------------------------------------------------------------------------
test_that("create massFlowTemplate class object and fill it with aligned samples", {
  tmp <- buildTMP(file = meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir)
  
  ## using default parameters
  tmp <- loadALIGNED(file = file.path(data_dir, "aligned.csv"),
                     template = file.path(data_dir, "template.csv"))
  expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == 2, tmp@params$bins == 0.05, tmp@params$cutoff == 0, tmp@params$realtime == FALSE))
  
  ## using selected parameters
  tmp <- loadALIGNED(file = file.path(data_dir, "aligned.csv"),
                     template = file.path(data_dir, "template.csv"),
                     rt_err = rt_err, realtime = TRUE, cutoff = cutoff)
  expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == rt_err, tmp@params$bins == 0.05, tmp@params$cutoff == cutoff, tmp@params$realtime == TRUE))
  
  expect_false(peaksVALIDATED(tmp))
  expect_true(validObject(tmp))
  expect_equal(rt_err, tmp@params$rt_err)
}) 


# buildANNO ---------------------------------------------------------------
test_that("create massFlowAnno class object ", {
  anno <- buildANNO(ds_file = ds_file, meta_file = meta_file, out_dir = anno_dir)
  ## Check object slots
  expect_true(class(anno) == "massFlowAnno")
  ## Check object slot values
  expect_true(anno@filepath == ds_file)
  expect_equal(anno@samples, read.csv(meta_file, stringsAsFactors = FALSE))
  expect_equal(anno@data, read.csv(ds_file, stringsAsFactors = FALSE))
  expect_true(class(anno@db) == "data.frame")
  expect_true(class(anno@anno) == "data.frame")
  expect_true(class(anno@mat) == "matrix")
  expect_true(class(anno@params) == "list")
})


