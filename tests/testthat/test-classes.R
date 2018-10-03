context("data classes")

test_that("initialise massFlowDB class object", {

  db <- buildDB(file = db_fname)
  db_table <- read.csv(file = db_fname, header = T, stringsAsFactors = F)

  # Check object values
  expect_true(class(db) == "massFlowDB")
  expect_equal(db@db, db_table)

  # Check setters
  expect_equal(chemicals(db), unique(db_table$dbname))

})

test_that("initialise massFlowTemplate class object", {

  tmp <- buildTMP(file = experiment_file)

  # Check object values
  expect_true(class(tmp) == "massFlowTemplate")
  expect_equal(tmp@samples$filepaths, experiment$filepaths)

})

test_that("initialise massFlowTemplate class object with DB", {

  db <- buildDB(file = db_fname)
  tmp <- buildTMP(file = experiment_file, db = db, rt_err = 10)

  # Check object values
  expect_true(class(tmp) == "massFlowTemplate")
  expect_equal(tmp@samples$filepaths, experiment$filepaths)

  expected_peakgrs <- data.frame(peakgr = c(1,8,10),
                                 tmp_peakgr = c(1,2,3),
                                 stringsAsFactors = F)
  obtained_peakgrs <- tmp@data[[1]] %>%
    filter(!is.na(cos)) %>%
    distinct(tmp_peakgr, peakgr) %>%
    arrange(tmp_peakgr)

  expect_equal(expected_peakgrs, obtained_peakgrs)
})


