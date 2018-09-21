context("data classes")

test_that("initialise massFlowDB class object", {

  db <- massFlowDB(file = db_fname)
  db_table <- read.csv(file = db_fname, header = T, stringsAsFactors = F)

  # Check object values
  expect_true(class(db) == "massFlowDB")
  expect_equal(db@db, db_table)

  # Check setters
  expect_equal(chemicals(db), unique(db_table$dbname))

})

test_that("initialise massFlowTemplate class object", {

  tmp <- massFlowTemplate(file = study_files)
  studyfiles <- read.csv(file = study_files, header = T, stringsAsFactors = F)
  tmp_table <- read.csv(studyfiles$filepaths[1], header = T, stringsAsFactors = F)

  # Check object values
  expect_true(class(tmp) == "massFlowTemplate")
  expect_equal(tmp@samples$filepaths, studyfiles$filepaths)

})

test_that("initialise massFlowTemplate class object with DB", {

  db <- massFlowDB(file = db_fname)
  tmp <- massFlowTemplate(file = study_files, db = db, rt_err = 10)
  studyfiles <- read.csv(file = study_files, header = T, stringsAsFactors = F)
  tmp_table <- read.csv(studyfiles$filepaths[1], header = T, stringsAsFactors = F)

  # Check object values
  expect_true(class(tmp) == "massFlowTemplate")
  expect_equal(tmp@samples$filepaths, studyfiles$filepaths)

  expected_peakgrs <- data.frame(peakgr = c(1,2,3),
                                 doi_peakgr = c(1,8,10), stringsAsFactors = F)
  obtained_peakgrs <- tmp@tmp %>%
    filter(!is.na(cos)) %>%
    distinct(peakgr, doi_peakgr) %>%
    arrange(peakgr)

  expect_equal(expected_peakgrs, obtained_peakgrs)
})


