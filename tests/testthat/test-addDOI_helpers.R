context("addDOI_helpers")

# checkFILE---------------------------------------------------------------------------------------------------
test_that("Checking grouped peaks csv tables via checkFILE is correct", {
  ## correct file
  checked_file <- massFlowR:::checkFILE(file = grouped_fnames[[1]])
  correct_file <-
    read.csv(file = grouped_fnames[[1]], stringsAsFactors = F)
  expect_equal(checked_file, correct_file)
  
  ## wrong filepath
  wrong_filepath <- gsub(".csv", "", grouped_fnames[[1]])
  expect_error(
    massFlowR:::checkFILE(file = wrong_filepath),
    paste("incorrect filepath for:", wrong_filepath)
  )
  
  ## wrong column names
  wrong_file <-
    setNames(correct_file, nm = c("peakID", "MZ", names(correct_file)[3:ncol(correct_file)]))
  wrong_file_filepath <-
    gsub(".csv", "wrong.csv", grouped_fnames[[1]])
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
  dt <- massFlowR:::checkFILE(file = grouped_fnames[[1]])
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
            ordered <-
              order(table(single_table$peakgr), decreasing = T)
            ordered <- factor(single_table$peakgr, levels = ordered)
            table_ordered <- single_table[order(ordered), ]
            
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

# addERRS ---------------------------------------------------------------------------------------------------------
test_that("addERR adds mz/rt windows correctly", {
  dt <- test_pks_rd
  dt[, c("mz_l",
         "mz_h",
         "rt_l",
         "rt_h")] <- c(dt$"mz" - mz_err,
                       dt$"mz" + mz_err,
                       dt$"rt" - rt_err,
                       dt$"rt" + rt_err)
  addERRS_out <-
    addERRS(dt = test_pks_rd,
            mz_err = mz_err,
            rt_err = rt_err)
  expect_identical(dt, addERRS_out)
  expect_true(addERRS_out[1, "mz_l"] == (test_pks_rd[1, "mz"] - mz_err))
})

# matchPEAK -------------------------------------------------------------------------------------------------------
test_that("matchPEAK", {
  ####---- match a table againts itself
  dt <- addERRS(dt = single_table,
                mz_err = mz_err,
                rt_err = rt_err)
  matchPEAK_out <- apply(dt, 1, FUN = matchPEAK, tmp = dt)
  matchPEAK_out <- do.call("rbindCLEAN", matchPEAK_out)
  
  expect_true(all(unlist(lapply(unique(single_table$peakid), function(peak) {
    peak_matched <-
      matchPEAK_out[matchPEAK_out$target_peakid == peak, "peakid"]
    peak_matched == peak
  }))))
  expect_true(all(unlist(lapply(unique(single_table$peakid), function(peak) {
    peak_matched <-
      matchPEAK_out[matchPEAK_out$target_peakid == peak, "mz"]
    peak_target <- single_table[single_table$peakid == peak, "mz"]
    peak_matched == peak_target
  }))))
  
  ####---- match against table with multiple matches
  dt1 <- addERRS(dt = single_table,
                 mz_err = mz_err,
                 rt_err = rt_err)
  peak_n <- 100
  dt2 <- bind_rows(dt1, dt1[peak_n,])
  
  matchPEAK_out <- apply(dt1, 1, FUN = matchPEAK, tmp = dt2)
  matchPEAK_out <- do.call("rbindCLEAN", matchPEAK_out)
  expect_true(matchPEAK_out$peakid[which(table(matchPEAK_out$peakid) > 1)] == peak_n)
  expect_true(table(matchPEAK_out$peakid)[peak_n] == 2)
  
  ####---- match against table without matches for selected peak
  dt1 <- addERRS(dt = single_table,
                 mz_err = mz_err,
                 rt_err = rt_err)
  peak_n <- 100
  dt2 <- dt1[dt1$peakid != peak_n,]
  matchPEAK_out <- apply(dt1, 1, FUN = matchPEAK, tmp = dt2)
  matchPEAK_out <- do.call("rbindCLEAN", matchPEAK_out)
  expect_false(peak_n %in% matchPEAK_out$target_peakid)
  expect_equal(nrow(matchPEAK_out), nrow(dt2))
})

# getCOS ----------------------------------------------------------------------------------------------------------
test_that("getCOS calculates cosines correctly", {

  ####---- identical peakgrs
  dt1 <- single_table
  dt2 <- single_table
  mat <- dt1 %>%
    mutate(target_peakgr = peakgr)
  target <- dt1
  target_peakgrs <- unique(dt2$peakgr)
  cos <- lapply(
    target_peakgrs,
    FUN = getCOS,
    target = target,
    mat = mat,
    tmp = dt1,
    bins = bins
  )
  cos <- do.call("rbindCLEAN", cos)
  expect_true(nrow(cos) == length(target_peakgrs))
  expect_true(all(round(cos$cos, digits = 10) == 1))
  
  ####---- compare the noisy peakgr and its original peakgr between two almost identical tables
  dt1 <-
    read.csv(noisy_fnames[1],
             header = T,
             stringsAsFactors = F)
  dt2 <-
    read.csv(noisy_fnames[2],
             header = T,
             stringsAsFactors = F)
  target <- dt2[dt2$peakgr == noisy_pkg,]
  mat <- dt1[dt1$peakgr == biggest_pkg,]
  mat$target_peakgr <- noisy_pkg
  getCOS_out <-
    getCOS(
      t_peakgr = noisy_pkg,
      target = target,
      mat = mat,
      tmp = dt1,
      bins = bins
    )
  expect_true(class(getCOS_out) == "data.frame")
  expect_true(getCOS_out$target_peakgr == noisy_pkg)
  expect_true(getCOS_out$peakgr == biggest_pkg)
  expect_true(getCOS_out$cos < 1)
  
  ####---- when no exact peak matches for the target peakgr are available
  ## this situation can happen if peakgrcls was selected with a peakgr that do not have exact matches to the target
  t_peakgr <- 1
  dt1 <-
    read.csv(noisy_fnames[1],
             header = T,
             stringsAsFactors = F)
  dt1 <- dt1[dt1$peakgr != t_peakgr,]
  dt2 <-
    read.csv(noisy_fnames[2],
             header = T,
             stringsAsFactors = F)
  target <- dt2[dt2$peakgr == t_peakgr,]
  ## create empty mat data frame which would normally contain peakgrs/peakgrcls peaks
  mat <- dt1[dt1$peakgr == t_peakgr,] %>%
    mutate(target_peakgr = t_peakgr)
  getCOS_out <-
    getCOS(
      t_peakgr = t_peakgr,
      target = target,
      mat = mat,
      tmp = dt1,
      bins = bins
    )
  expect_true(is.null(getCOS_out))
  
  ####---- two different tables
  dt1 <-
    read.csv(grouped_fnames[1],
             header = T,
             stringsAsFactors = F)
  dt2 <-
    read.csv(grouped_fnames[2],
             header = T,
             stringsAsFactors = F)
  dt1 <- addERRS(dt = dt1,
                 mz_err = mz_err,
                 rt_err = rt_err)
  dt2 <- addERRS(dt = dt2,
                 mz_err = mz_err,
                 rt_err = rt_err)
  t_peakgr <- 26
  target <- dt2[dt2$peakgr == t_peakgr,]
  mat <- apply(target, 1, FUN = matchPEAK, tmp = dt1)
  mat <- do.call("rbindCLEAN", mat)
  getCOS_out <-
    getCOS(
      t_peakgr = t_peakgr,
      target = target,
      mat = mat,
      tmp = dt1,
      bins = bins
    )
  expect_true(class(getCOS_out) == "data.frame")
  expect_true(round(getCOS_out$cos, digits = 6) == 0.999552)
  
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
