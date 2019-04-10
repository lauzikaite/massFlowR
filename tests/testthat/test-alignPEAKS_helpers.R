context("helper functions for alignPEAKS method")

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

# getRTbins ---------------------------------------------------------------------------------------------------------
test_that("getRTbins() correctly splits two dataframes to rt bins", {
  ## make dummy table:
  ## 20 peaks, 4 peakgroups: peakgroups number 2 and 3 are very close
  ## split dummy table into two rt bins
  ## as number 2 and 3 are very close, they will be included in both bin 1 and 2
  dummy <- data.frame(peakid = 1:20,
                      peakgr = rep(1:4, each = 5),
                      rt = rep(c(100, 299, 301, 400), each = 5))
  ncores <- 2
  ## use the same table as both ds and tmp
  out <-
    getRTbins(
      ds = dummy,
      tmp = dummy,
      ds_var_name = "peakgr",
      tmp_var_name = "peakgr",
      mz_err = mz_err,
      rt_err = rt_err,
      ncores = ncores
    )

  ## number of created rt bins equals ncores
  expect_true(all(length(out$ds) == ncores, length(out$tmp) == ncores))
  
  ## all peakids are exported
  expect_true(all(dummy$peakid %in% c(out$tmp[[1]]$peakid, out$tmp[[2]]$peakid)))
  expect_true(all(dummy$peakid %in% c(out$ds[[1]]$peakid, out$ds[[2]]$peakid)))
  
  ## since ds and tmp are identical, all peaks in the ds bin must also be present in the corresponding tmp bin
  expect_true(all(
    out$ds[[1]]$peakid %in% out$tmp[[1]]$peakid,
    out$ds[[2]]$peakid %in% out$tmp[[2]]$peakid
  ))
  
  ## obtained tmp bins must overlap for common rt values
  rt_bin_1 <- c(min(out$tmp[[1]]$rt - rt_err),
                max(out$tmp[[1]]$rt + rt_err))
  rt_bin_2 <- c(min(out$tmp[[2]]$rt - rt_err),
                max(out$tmp[[2]]$rt + rt_err))
  peakgrs_in_common <-
    dummy$peakgr[which((dummy$rt - rt_err) >= rt_bin_2[1] &
                              (dummy$rt + rt_err) <= rt_bin_1[2])]
  peaks_in_common <- dummy$peakid[which(dummy$peakgr %in% peakgrs_in_common)]
  expect_identical(peaks_in_common,
                   base::intersect(out$tmp[[1]]$peakid, out$tmp[[2]]$peakid))
})

# orderBYrt ---------------------------------------------------------------------------------------------------------
test_that("orderBYrt() correctly orders dataset using rt", {
  ## create sample dataset
  dt <- data.frame(
    peakid = c(3, 9, 21, 36, 46, 78, 81, 93, 124, 149, 4, 17, 52, 62, 67, 119, 123, 6, 19, 20, 60, 101, 186, 13, 106, 117, 118, 1, 10, 35, 2, 12, 50, 5, 18, 63, 8, 22, 68, 11, 33, 87, 7, 24),
    mz = c(343, 365, 344, 366, 381, 345, 388.100006103516, 362.100006103516, 382, 444.100006103516, 524.200012207031, 525.200012207031, 526.200012207031, 249.100006103516, 255.199996948242, 282.200012207031, 309.200012207031, 522.200012207031, 531.200012207031, 523.200012207031, 389.200012207031, 344.200012207031, 557.100036621094, 360, 431.100006103516, 370.100006103516, 416.100006103516, 508.200012207031, 509.200012207031, 367.200012207031, 496.200012207031, 497.200012207031, 498.200012207031, 526.100036621094, 527.100036621094, 385.100006103516, 522.200012207031, 523.200012207031, 524.200012207031, 496.200012207031, 497.200012207031, 498.200012207031, 502.100006103516, 503.100006103516),
    rt = c(2678.218, 2679.783, 2678.218, 2679.783, 2679.783, 2678.218, 2678.218, 2678.218, 2676.653, 2679.783, 3662.573, 3662.573, 3662.573, 3662.573, 3661.009, 3662.573, 3661.009, 3344.888, 3344.888, 3344.888, 3343.323, 3344.888, 3343.323, 2684.478, 2684.478, 2682.913, 2684.478, 3515.468, 3517.033, 3517.033, 3384.012, 3384.012, 3384.012, 3168.048, 3168.048, 3168.048, 3409.051, 3409.051, 3409.051, 3316.719, 3316.719, 3316.719, 3157.093, 3157.093),
    into = c(26239739.4908571, 15389161.888963, 5872251.05168, 3427455.12957693, 2193772.47438462, 945843.066319147, 866478.105453334, 733169.002746268, 463242.195000001, 277061.806153847, 26223613.1127728, 7829455.07127275, 1933029.53740001, 1307022.7299762, 1168879.00468966, 515198.332384615, 492462.693627908, 23482537.5170256, 7268277.8061628, 6909134.43804878, 1374905.10836585, 672994.049432435, 122603.611878788, 10062819.6947222, 637630.751882354, 525470.723055555, 523740.547333333, 53910492.9448513, 15123847.3582767, 3493079.33254168, 38200390.2189089, 10219380.1002272, 1972152.95488889, 25760336.3028836, 7615134.77321737, 1269534.71757576, 19327938.1211429, 5775851.12990477, 1156396.91788889, 13639592.8685714, 3764370.88609524, 780206.794606061, 20113368.9417143, 5122895.84110345),
    peakgr = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 1, 1, 1, 2, 2, 2, 5, 5, 5, 8, 8, 8, 9, 9, 9, 7, 7)
    )
  ## get peak-group order by their median rt
  dt_vars <- unique(dt$peakgr)
  dt_rt_med <- sapply(dt_vars, function(x) {
    median(dt[dt$peakgr == x, "rt"])
  })
  expect_true(all(sapply(2:length(dt_rt_med), function(x) {
    dt_rt_med[order(dt_rt_med)][x] >  dt_rt_med[order(dt_rt_med)][x-1]
  })))
  dt_var_levels <- factor(dt[, "peakgr"], levels = dt_vars[order(dt_rt_med)])
  dt_expected <- dt[order(dt_var_levels),]
  
  ## run function
  dt_ordered <- orderBYrt(dt = dt, var_name = "peakgr")
  
  ## compare
  expect_identical(dt_expected, dt_ordered)
  dt_rt_ordered_med <- sapply(unique(dt_ordered$peakgr), function(x) {
    median(dt_ordered[dt_ordered$peakgr == x, "rt"])
  })
  expect_true(all(sapply(2:length(dt_rt_ordered_med), function(x) {
    dt_rt_ordered_med[x] > dt_rt_ordered_med[x - 1]
  })))
  expect_true(all(sapply(unique(dt_ordered$peakgr), function(x) {
    all(diff(order(dt_ordered[dt_ordered$peakgr == x, "peakid"])) == 1)
  })))
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

# getCOSmat ---------------------------------------------------------------------------------------------------------
test_that("getCOSmat() returns correct cosine matrix", {
  
  ####---- compare the noisy peakgr and its original peakgr between two ALMOST identical tables ----####
  tmp <- read.csv(noisy_fnames[1],
                  header = T,
                  stringsAsFactors = F)
  doi <- read.csv(noisy_fnames[2], 
                  header = T,
                  stringsAsFactors = F)
  rt_bins <- getRTbins(ds = doi,
                       tmp = tmp,
                       ds_var_name = "peakgr",
                       tmp_var_name = "peakgr",
                       mz_err = mz_err,
                       rt_err = rt_err,
                       ncores = 2)
  doi_bins <- rt_bins$ds
  tmp_bins <- rt_bins$tmp
  doi_vars_by_bins <- lapply(doi_bins, function(x) unique(x$peakgr))
  bin <- which(sapply(doi_vars_by_bins, function(x) noisy_pkg %in% x))
  ds_bin <- doi_bins[[bin]]
  tmp_bin <- tmp_bins[[bin]]
  cos_matches <- getCOSmat(
    ds_bin = doi_bins[[bin]],
    ds_var = "peakgr",
    tmp_bin = tmp_bins[[bin]],
    tmp_var = "peakgr",
    mz_err = mz_err,
    rt_err = rt_err,
    bins = bins
  )
  cos_mat <- cos_matches[[1]]
  ####---- 30 peak-group pairs have cos of 1 
  expect_true(table(cos_mat[which(cos_mat > 0)])["1"] == 30)
  
  ####---- those that have cos < 1 are the noisy peakgroups
  r_names <- names(Filter(length, apply(cos_mat, 1, function(x) {
    which(x > 0 & x < 1)
  })))
  c_names <- names(Filter(length, apply(cos_mat, 2, function(x) {
    which(x > 0 & x < 1)
  })))
  
  ## in tmp, the peakgroup from which two noisy peakgrs were generated in the later file is biggest_pkg (ie 19)
  ## in doi, the noisy peakgroups generated from tmp peakgr are noisy_pkg (ie 63) and biggest_pkg( ie 19)
  ## cos mat rownames correspond to tmp peakgroups and colnames correspond to doi peakgroups
  expect_true(r_names == biggest_pkg)
  expect_true(all(c_names == c(noisy_pkg, biggest_pkg)))
  
  ####---- identical tables ----####
  doi <- read.csv(grouped_fnames[1],
                  header = T,
                  stringsAsFactors = F)
  tmp <- read.csv(grouped_fnames[1],
                  header = T,
                  stringsAsFactors = F)
  tmp <- tmp[, c("peakid",	"mz", "rt", "into", "peakgr")]
  rt_bins <- getRTbins(ds = doi,
                       tmp = tmp,
                       ds_var_name = "peakgr",
                       tmp_var_name = "peakgr",
                       mz_err = mz_err,
                       rt_err = rt_err,
                       ncores = 1)
  doi_bins <- rt_bins$ds
  tmp_bins <- rt_bins$tmp
  doi_vars_by_bins <- lapply(doi_bins, function(x) unique(x$peakgr))
  cos_matches <- getCOSmat(
    ds_bin = doi_bins[[bin]],
    ds_var = "peakgr",
    tmp_bin = tmp_bins[[bin]],
    tmp_var = "peakgr",
    mz_err = mz_err,
    rt_err = rt_err,
    bins = bins
  )
  cos_mat <- cos_matches[[1]]
  ## since doi and tmp are identical, all cos for corresponding peakgroups should be 1
  expect_true(table(cos_mat[which(cos_mat > 0)])["1"] == length(unique(doi$peakgr)))
  r_names <- names(Filter(length, apply(cos_mat, 1, function(x) {
    which(x > 0)
  })))
  c_names <- names(Filter(length, apply(cos_mat, 2, function(x) {
    which(x > 0)
  })))
  expect_true(all(r_names == c_names))
  
  ####---- different tables ----####
  tmp <- read.csv(grouped_fnames[1],
                  header = T,
                  stringsAsFactors = F)
  doi <- read.csv(grouped_fnames[2],
                  header = T,
                  stringsAsFactors = F)
  tmp <- tmp[, c("peakid",	"mz", "rt", "into", "peakgr")]
  ncores <- 1
  rt_bins <- getRTbins(ds = doi,
                       tmp = tmp,
                       ds_var_name = "peakgr",
                       tmp_var_name = "peakgr",
                       mz_err = mz_err,
                       rt_err = rt_err,
                       ncores = ncores)
  doi_bins <- rt_bins$ds
  tmp_bins <- rt_bins$tmp
  doi_vars_by_bins <- lapply(doi_bins, function(x) unique(x$peakgr))
  cos_matches <- list(
    getCOSmat(
      ds_bin = doi_bins[[ncores]],
      ds_var = "peakgr",
      tmp_bin = tmp_bins[[ncores]],
      tmp_var = "peakgr",
      mz_err = mz_err,
      rt_err = rt_err,
      bins = bins)
  )
  expect_equal(length(cos_matches[[ncores]][[2]]), length(unique(doi$peakgr)))
  expect_equal(nrow(cos_matches[[ncores]][[1]]), length(unique(tmp$peakgr)))
  expect_equal(ncol(cos_matches[[ncores]][[1]]), length(unique(doi$peakgr)))
  expect_equal(max(cos_matches[[ncores]][[1]]), 0.9999)
  expect_true(cos_matches[[ncores]][[1]]["22","27"] == max(cos_matches[[ncores]][[1]]))
})

# getMATCHES ---------------------------------------------------------------------------------------------------------
test_that("getMATCHES() ", {
  tmp <- addERRS(single_table, mz_err = mz_err, rt_err = rt_err)
  ## target peaks are: 
  ## 1) a peak from the same table
  ## 2) a modified version of the same peak, within matching range 
  ## 3) a modified version of the same peak, outside of matching range 
  target_peak1 <- tmp[1, ]
  target_peak2 <- tmp[1, ]
  target_peak2$mz <- target_peak2$mz - (mz_err/2)
  target_peak2$rt <- target_peak2$rt + (rt_err/2)
  target_peak3 <- tmp[1,]
  target_peak3$mz <- target_peak3$mz - (mz_err*3)
  target_peak3$rt <- target_peak3$rt + (rt_err*3)
  
  target_peaks <- rbind(target_peak1, target_peak2, target_peak3)
  target_peaks <- addERRS(target_peaks, mz_err = mz_err, rt_err = rt_err)
  matches <- lapply(1:nrow(target_peaks), getMATCHES, target = target_peaks, tmp = tmp, tmp_var = "peakgr", target_var = "peakgr")
  
  expect_true(length(matches) == nrow(target_peaks))
  expect_equal(sapply(matches, nrow), list(1, 1, NULL))
  expect_equal(sapply(matches, "[[", "peakid"), list(tmp[1, "peakid"], tmp[1, "peakid"], NULL))
  expect_equal(sapply(matches, colnames),
               c(rep(list(
                 c(
                   "peakid",
                   "mz",
                   "rt",
                   "into",
                   "tmp_var",
                   "target_peakid",
                   "target_var"
                 )
               ), 2), list(NULL)))
})

# buildVECTOR ---------------------------------------------------------------------------------------------------------
test_that("buildVECTOR()", {
  
  spec <- data.frame(mz = seq(from = 385.1000061,
                              to = 527.1000366,
                              by = bins))
  spec$bin <- 1:nrow(spec)
  peaks <- data.frame(
    mz = c(526.1000366,
           527.1000366,
           386.1000061),
    into = c(27029849.76,
             7613856.402,
             361034.734)
  )
  sample_vec <- buildVECTOR(spec = spec, peaks = peaks)
  expect_true(all(which(sample_vec > 0) == c(101, 14101, 14201)))
  sum(sample_vec[which(sample_vec > 0)])
})

# scaleSPEC ---------------------------------------------------------------------------------------------------------
test_that("scaleSPEC() ", {
  ## Version A - scale to unit length
  spec <- data.frame(into = c(
    rep(0, 10),
    0.0128555075994933,
    rep(0, 10),
    0.271109618049552,
    rep(0, 10),
    0.9624626283266
  ))
  spec_scaled <- scaleSPEC(spec)
  expect_identical(length(spec_scaled), nrow(spec))
  expect_identical(length(which(spec_scaled > 0)), length(which(spec$into > 0)))
  spec_length <- sum(spec_scaled * spec_scaled)
  expect_identical(spec_length, 1)
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