context("data classes and classes constructors")

test_that("initialise massFlowTemplate class object using default parameters", {

  tmp <- buildTMP(file = meta_fname, out_dir = data_dir)
  data1 <- tmp@data[[1]]
  
  # Check object values
  expect_true(class(tmp) == "massFlowTemplate")
  expect_equal(tmp@samples$filepaths, metadata$filepaths)
  expect_equal(tmp@samples[which(tmp@samples$aligned == TRUE),"filepaths"], metadata$filepaths[1])
  expect_true(all(data1$peakid %in% single_table$peakid))
  expect_true(all_equal(
    data1[match(single_table$peakid, data1$peakid),c(names(single_table))],
    single_table))
  expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == 2, tmp@params$bins == 0.01))
  
  # Check getters
  expect_true(filepath(tmp) == meta_fname)

})

test_that("create massFlowTemplate class object and fill it with aligned samples", {

  tmp <- buildTMP(file = meta_fname, out_dir = data_dir, rt_err = rt_err)
  tmp <- alignPEAKS(tmp, out_dir = data_dir, write_int = F)
  tmp <- loadALIGNED(file = file.path(data_dir, "aligned.csv"),
                     template = file.path(data_dir, "template.csv"),
                     rt_err = rt_err)
  
  expect_false(peaksVALIDATED(tmp))
  expect_true(validObject(tmp))
  expect_equal(rt_err, tmp@params$rt_err)
  
}) 

