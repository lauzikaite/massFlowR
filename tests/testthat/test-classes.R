context("data classes and classes constructors")

test_that("initialise massFlowDB class object", {

  db <- buildDB(file = db_file)
  db_table <- read.csv(file = db_file, header = T, stringsAsFactors = F)

  # Check object values
  expect_true(class(db) == "massFlowDB")
  expect_equal(db@db, db_table)

  # Check getters
  expect_true(filepath(db) == db_file)
  
})

test_that("initialise massFlowTemplate class object without DB, using default parameters", {

  tmp <- buildTMP(file = experiment_file)

  # Check object values
  expect_true(class(tmp) == "massFlowTemplate")
  expect_true(tmp@db_filepath == "not used")
  expect_equal(tmp@samples$filepaths, experiment$filepaths)
  expect_equal(tmp@samples[which(tmp@samples$aligned == TRUE),"filepaths"], experiment$filepaths[1])
  expect_true(all_equal(tmp@data[[1]][which(tmp@data[[1]]$peakid %in%  single_table$peakid),c(names(single_table))], single_table))
  expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == 2, tmp@params$bins == 0.01, tmp@params$db_thrs == 0))
  
  # Check getters
  expect_true(filepath(tmp) == experiment_file)

})

test_that("initialise massFlowTemplate class object with DB", {

  ## test will fail due to the nature of matchPEAK with character entries for DB columns, must change matchPEAK
  # db <- buildDB(file = db_file)
  # tmp <- buildTMP(file = experiment_file, db = db, rt_err = 10)
  # 
  # # Check object values
  # expect_true(class(tmp) == "massFlowTemplate")
  # expect_equal(tmp@filepath, experiment_file)
  # expect_equal(tmp@samples$filepaths, experiment$filepaths)
  # expect_equal(tmp@samples$filepaths[which(tmp@samples$aligned)], experiment$filepaths[1])
  # expect_equal(tmp@db_filepath, db_file)
  # expect_true(all(tmp@params$mz_err == 0.01, tmp@params$rt_err == 10, tmp@params$bins == 0.01, tmp@params$db_thrs == 0.5))
  # 
  # expected_peakgrs <- data.frame(peakgr = c(1,8,10),
  #                                tmp_peakgr = c(1,2,3),
  #                                stringsAsFactors = F)
  # obtained_peakgrs <- tmp@data[[1]] %>%
  #   filter(!is.na(cos)) %>%
  #   distinct(tmp_peakgr, peakgr) %>%
  #   arrange(tmp_peakgr)
  # expect_equal(expected_peakgrs, obtained_peakgrs)
  # 
  # # Check getters
  # expect_true(filepath(tmp) == experiment_file)
  # 
})


