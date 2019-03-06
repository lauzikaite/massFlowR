# rankCOS -------------------------------------------------------------------------------------------------
test_that("", {
 
  ## ver1
  cos_list <- c(0, 0,
                0.9, 0.8,
                0.85, 0.75)
  cos <- matrix(cos_list, nrow = 2, ncol = 3)
  
  ## ver2
  cos_list <- c(0, 0, 0, 0,
                0, 0.9, 0.8, 0.5,
                0.5, 0, 0, 0.4,
                0, 0.8, 0.7, 0.6)
  cos <- matrix(cos_list, nrow = 4, ncol = 4)
  
  ## version4, column 7 will get assigned with its 2nd best option - row 3 (for which it is the first best option)
  cos_list <- c(0, 0, 0, 0,
                0, 0.9, 0.7, 0,
                0, 0.1, 0, 0.9,
                0, 0.1, 0.6, 0.8,
                0.5, 0, 0.5, 0,
                0.1, 0.2, 0, 0.4,
                0, 0.85, 0.8, 0)
  cos <- matrix(cos_list, nrow = 4, ncol = 7)
  
  ## version4, column 4 will get assigned with its second best option - row 3 (for which it is the third best option)
  cos_list <- c(0, 0, 0, 0, 0, 
                0, 0.9, 0.7, 0, 0,
                0, 0.1, 0, 0.9, 0,
                0, 0.1, 0.6, 0.8, 0.5,
                0.5, 0, 0.5, 0, 0,
                0.1, 0.2, 0, 0.4, 0, 
                0, 0.85, 0.8, 0, 0.9)
  cos <- matrix(cos_list, nrow = 5, ncol = 7)
  

  
})