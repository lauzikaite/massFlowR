context("main groupPEAKS function")

# groupPEAKS --------------------------------------------------------------
test_that("groupPEAKS() correctly recognise wrong inputs", {
  cwt_test <- cwt
  cwt_test@mzdiff <- -0.05
  
  expect_error(groupPEAKS(file = "not/a/file", out_dir = data_dir),
               regexp = "incorrect filepath for 'file' provided")
  expect_error(groupPEAKS(file = meta_fname, out_dir = "Not/A/Path"),
                 regexp = "incorrect filepath for 'out_dir' provided")
  expect_error(groupPEAKS(file = meta_fname, out_dir = data_dir),
               regexp = "'cwt' has to be a 'CentWaveParam' object")
  expect_error(groupPEAKS(file = grouped_fnames[1], out_dir = data_dir, cwt = cwt_test),
               regexp = "'files' table must contain columns: filename, run_order, raw_filepath")
  expect_warning(groupPEAKS(file = meta_fname, out_dir = data_dir, cwt = cwt_test),
                 regexp =  "mzdiff of -0.05 was selected. Switching mzdiff to 0 ...")
})
  
