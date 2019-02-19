context("data classes and classes constructors")

test_that("initialise massFlowDB class object", {

  db <- buildDB(file = db_fname)
  db_table <- read.csv(file = db_fname, header = T, stringsAsFactors = F)

  # Check object values
  expect_true(class(db) == "massFlowDB")
  expect_equal(db@db, db_table)

  # Check getters
  expect_true(filepath(db) == db_fname)
  
})

test_that("initialise massFlowTemplate class object using default parameters", {

  tmp <- buildTMP(file = meta_fname, out_dir = data_dir)

  # Check object values
  expect_true(class(tmp) == "massFlowTemplate")
  expect_equal(tmp@samples$filepaths, metadata$filepaths)
  expect_equal(tmp@samples[which(tmp@samples$aligned == TRUE),"filepaths"], metadata$filepaths[1])
  expect_true(all_equal(tmp@data[[1]][which(tmp@data[[1]]$peakid %in%  single_table$peakid),c(names(single_table))], single_table))
  expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == 2, tmp@params$bins == 0.01))
  
  # Check getters
  expect_true(filepath(tmp) == meta_fname)

})

test_that("create massFlowTemplate class object and fill it with aligned samples", {
  
  
}) 

