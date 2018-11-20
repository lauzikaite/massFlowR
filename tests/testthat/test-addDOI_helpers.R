context("addDOI_helpers")

# checkFILE---------------------------------------------------------------------------------------------------
test_that("Checking grouped peaks csv tables via checkFILE is correct", {
  ## correct file
  checked_file <- massFlowR:::checkFILE(file = grouped_files[[1]])
  correct_file <-
    read.csv(file = grouped_files[[1]], stringsAsFactors = F)
  expect_equal(checked_file, correct_file)
  
  ## wrong filepath
  wrong_filepath <- gsub(".csv", "", grouped_files[[1]])
  expect_error(
    massFlowR:::checkFILE(file = wrong_filepath),
    paste("incorrect filepath for:", wrong_filepath)
  )
  
  ## wrong column names
  wrong_file <-
    setNames(correct_file, nm = c("peakID", "MZ", names(correct_file)[3:ncol(correct_file)]))
  wrong_file_filepath <-
    gsub(".csv", "wrong.csv", grouped_files[[1]])
  write.csv(wrong_file, file = wrong_file_filepath, row.names = F)
  expected_error <-
    paste("incorrect file:",
          wrong_file_filepath,
          "\n",
          "missing columns: peakid, mz")
  expect_error(massFlowR:::checkFILE(file = wrong_file_filepath),
               expected_error)
  unlink(x = wrong_file_filepath)
})


# getCLUSTS -------------------------------------------------------------------------------------------------------
test_that("Clustering of peak-groups via getCLUSTS is correct", {
  dt <- massFlowR:::checkFILE(file = grouped_files[[1]])
  clustered_dt <- massFlowR:::getCLUSTS(dt = dt)
  
  expect_equal(nrow(clustered_dt), nrow(dt))
  expect_true(all(dt$peakid %in% clustered_dt$peakid))
  expect_true(all(dt$peakgr %in% clustered_dt$peakgr))
  ## check order of peak-groups
  expect_true(all(order(
    table(clustered_dt$peakgr)[unique(clustered_dt$peakgr)], decreasing = T
  ) == 1:max(dt$peakgr)))
  
})


# orderPEAKS ------------------------------------------------------------------------------------------------------
test_that("Peak-group ordering by complexity and intensity via orderPEAKS is correct",
          {
            ## according to how many peaks per peak-group
            ordered <- order(table(single_table$peakgr), decreasing = T)
            ordered <- factor(single_table$peakgr, levels = ordered)
            table_ordered <- single_table[order(ordered),]
            
            ## according to peak intensity (aka peakid)
            peaks_ordered <-
              sapply(unique(table_ordered$peakgr), function(pkg) {
                ord <-
                  order(table_ordered[which(table_ordered$peakgr == pkg), "peakid"])
                if (all(ord == 1:length(ord))) {
                  TRUE
                } else {
                  FALSE
                }
              })
            expect_true(all(unlist(peaks_ordered)))
            
            ## compare peakid order in orderPEAKS output
            orderPEAKS_out <- orderPEAKS(dt = single_table)
            expect_true(all(orderPEAKS_out$peakid == table_ordered$peakid))
            
          })


# matchPEAK -------------------------------------------------------------------------------------------------------
test_that("matchPEAK", {
  
  
  
  
})
  








# compareCOS -------------------------------------------------------------------------------------------------
test_that("Selection of top matches via compareCOS is correct", {
  ## sample table consists of 3 target peakgrs, each matched by the same 3 template peakgrs
  cos_list <- c(seq(
    from = 0.1,
    to = 1,
    length.out = 8
  ), 1)
  cos <- data.frame(
    target_peakgr = c(rep(1, 3), rep(2, 3), rep(3, 3)),
    peakgr = c(rep(c(1, 2, 3), 3)),
    cos = cos_list,
    stringsAsFactors = F
  )
  expect_error(top <-
                 massFlowR:::compareCOS(cos = cos),
               "identical cosines were found!")
  
  ## 1-1 is top pair by rank
  ## 2-2 is second best for target peakgr (choosing from what is left unassigned)
  ## 3-3 is third best for target peakgr (choosing from what is left unassigned)
  cos_list <- seq(from = 0.1,
                  to = 1,
                  length.out = 9)
  cos_list <- cos_list[c(9, 6, 3, 8, 5, 2, 7, 4, 1)]
  cos <- cos %>%
    mutate(cos = cos_list)
  
  top <- massFlowR:::compareCOS(cos = cos)
  expected_top <- data.frame(
    target_peakgr = c(1, 2, 3),
    peakgr = c(1, 2, 3),
    cos = cos_list[c(1, 5, 9)],
    top = rep(T, 3)
  )
  expect_equal(top, expected_top)
  
  ## when one of the target peakgrs doesn't have a cos >0
  cos_list <- c(rep(0, 3), seq(
    from = 0.1,
    to = 1,
    length.out = 6
  ))
  cos <- data.frame(
    target_peakgr = c(rep(1, 3), rep(2, 3), rep(3, 3)),
    peakgr = c(rep(c(1, 2, 3), 3)),
    cos = cos_list,
    stringsAsFactors = F
  )
  top <- massFlowR:::compareCOS(cos = cos)
  ## returned table will be ordered by target_peakgrs and their highest cosine
  expected_top <- data.frame(
    target_peakgr = c(3, 2),
    peakgr = c(3, 2),
    cos = c(cos_list[c(9, 5)]),
    top = c(T, T)
  )
  expect_equal(top, expected_top)
})
