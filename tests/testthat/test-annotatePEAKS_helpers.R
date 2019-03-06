context("annotatePEAKS_helpers")


# rankCOS ---------------------------------------------------------------------------------------------------------
test_that("rankCOS() correctly ranks cosines, giving 1 to the highest cosine", {
  cos_list <- c(0, 0, 0.9, 0,
                0.6, 0.9, 0.7, 0,
                0.7, 0, 0.8 , 0.6, 
                0.8, 0.85, 0, 0)
  cos <- matrix(cos_list, nrow = 4, ncol = 4)
  ## rank for columns
  cos_ranked <- apply(cos, 2, FUN = massFlowR:::rankCOS)
  expected <- matrix(c(0, 0, 1, 0,
                       3, 1, 2, 0,
                       2, 0, 1, 3,
                       2, 1, 0, 0),
                     nrow = 4, ncol = 4)
  expect_equal(cos_ranked, expected)
  
  ## rank for rows
  cos_ranked <- apply(cos, 1, FUN = massFlowR:::rankCOS)
  expected <- matrix(c(0, 3, 2, 1,
                       0, 1, 0, 2,
                       1, 3, 2, 0,
                       0, 0, 1, 0),
                     nrow = 4, ncol = 4)
  expect_equal(cos_ranked, expected)
  
})
  







# assignCOS -------------------------------------------------------------------------------------------------
test_that("assignCOS() correctly assigns pairs to maximise top-top matching", {
 
  ## version 1
  cos_list <- c(0, 0, 0, 0,
                0, 0.9, 0.8, 0.5,
                0.5, 0, 0, 0.4,
                0, 0.8, 0.7, 0.6)
  cos <- matrix(cos_list, nrow = 4, ncol = 4)
  cos_assigned <- massFlowR:::assignCOS(cos)
  expected <- matrix(c(rep(NA, 4),
                       NA, TRUE, NA, NA,
                       TRUE, rep(NA, 3),
                       rep(NA, 2), TRUE, NA),
                     nrow = 4, ncol = 4)
  expect_equal(cos_assigned, expected)
  
  ## version2, column 7 will get assigned with its 2nd best option - row 3 (for which it is the first best option)
  cos_list <- c(0, 0, 0, 0,
                0, 0.9, 0.7, 0,
                0, 0.1, 0, 0.9,
                0, 0.1, 0.6, 0.8,
                0.5, 0, 0.5, 0,
                0.1, 0.2, 0, 0.4,
                0, 0.85, 0.8, 0)
  cos <- matrix(cos_list, nrow = 4, ncol = 7)
  cos_assigned <- massFlowR:::assignCOS(cos)
  expected <- matrix(c(rep(NA, 4),
                       NA, TRUE, rep(NA, 2),
                       rep(NA, 3), TRUE,
                       rep(NA, 4),
                       TRUE, rep(NA, 3),
                       rep(NA, 4),
                       rep(NA, 2), TRUE, NA),
                     nrow = 4, ncol = 7)
  expect_equal(cos_assigned, expected)
  
  ## version3, column 4 will get assigned with its second best option - row 3 (for which it is the third best option)
  cos_list <- c(0, 0, 0, 0, 0, 
                0, 0.9, 0.7, 0, 0,
                0, 0.1, 0, 0.9, 0,
                0, 0.1, 0.6, 0.8, 0.5,
                0.5, 0, 0.5, 0, 0,
                0.1, 0.2, 0, 0.4, 0, 
                0, 0.85, 0.8, 0, 0.9)
  cos <- matrix(cos_list, nrow = 5, ncol = 7)
  cos_assigned <- massFlowR:::assignCOS(cos)
  expected <- matrix(c(rep(NA, 5),
                       NA, TRUE, NA, NA, NA,
                       rep(NA, 3), TRUE, NA,
                       rep(NA, 2), TRUE, rep(NA, 2),
                       TRUE, rep(NA, 4),
                       rep(NA, 5),
                       rep(NA, 4), TRUE),
                     nrow = 5, ncol = 7)
  expect_equal(cos_assigned, expected)
  
  ## version4, column 3 is assigned with its third best option, column 4 is assigned with its second best option
  cos_list <- c(0, 0, 0.9, 0,
                0.6, 0.9, 0.7, 0,
                0.7, 0, 0.8 , 0.6, 
                0.8, 0.85, 0, 0)
  cos <- matrix(cos_list, nrow = 4, ncol = 4)
  cos_assigned <- massFlowR:::assignCOS(cos)
  expected <- matrix(c(NA, NA, TRUE, NA,
                       NA, TRUE, NA, NA,
                       NA, NA, NA, TRUE,
                       TRUE, NA, NA, NA),
                     nrow = 4, ncol = 4)
  expect_equal(cos_assigned, expected)

})