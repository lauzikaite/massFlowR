context("helper functions for validPEAKS method")

####---- obtain object with different samples
tmp <- buildTMP(file = meta_fname, out_dir = data_dir, rt_err = rt_err)
tmp <- alignPEAKS(tmp, out_dir = data_dir)

# extractPEAKGR ---------------------------------------------------------------------------------------------------
test_that("extractPEAKGR extract intensities for each peakgr from each sample using massFlowTemplate object", {
  ## take a peakgr
  pkg <- 1
  extractPEAKGR_out <- extractPEAKGR(pkg = pkg, object = tmp)
  expect_equal(names(extractPEAKGR_out), c("peakid", "into", "maxo", "mz", "rt", "run_order", "peakgr"))
  expect_true(all(unique(extractPEAKGR_out$run_order) == 1:2))
  ## second sample has the same peaks as first sample + 1 unique for this peakgr
  expect_equal(nrow(extractPEAKGR_out), nrow(single_table[single_table$peakgr == pkg, ]) * 2 + 1)
  for (x in unique(extractPEAKGR_out$run_order)) {
    out_x <- tmp@data[[x]]
    expect_equal(as.list(out_x[out_x$peakgr == pkg, c("into", "maxo", "mz", "rt")]),
                 as.list(extractPEAKGR_out[extractPEAKGR_out$run_order == x, c("into", "maxo", "mz", "rt")]))
  }
  ## peakids in the extractPEAKGR_out are tmp peakids (not original sample peakids)
  out_x <- tmp@data[[2]]
  expect_equal(tmp@tmp[tmp@tmp$peakgr == pkg, "peakid"],
               extractPEAKGR_out[extractPEAKGR_out$run_order == x, c("peakid")])
})

# validPEAKGR ---------------------------------------------------------------------------------------------------
test_that("validPEAKGR returns correctly communities of peaks", {
  ####---- use object with different samples
  pkg <- 10 # take a peakgr from final template
  tmp_pkg <- tmp@tmp[tmp@tmp$peakgr == pkg, ]
  ## extract intensities 
  pkg_ints <- extractPEAKGR(pkg = pkg, object = tmp)
  ## obtain communities
  validPEAKGR_out <- validPEAKGR(pkg = pkg, pkg_ints = pkg_ints, out_dir = out_dir, min_samples_n = 0, cor_thr = 0.75)
  ## all peaks are listed in the exported communities
  expect_true(all(unlist(sapply(validPEAKGR_out, "[[", "peakid")) %in% tmp_pkg[ ,"peakid"]))
  ## exported into/mz/rt are correct and match up the original values for each sample
  ## for each community
  for (x in 1:length(validPEAKGR_out)) {
    out_x <- validPEAKGR_out[[x]]
    ## for each sample
    for (r in unique(out_x$run_order)) {
      out_x_r <- out_x[out_x$run_order == r, ]
      data_r <- tmp@data[[r]]
      expect_equal(
        as.list(out_x_r[match(data_r$tmp_peakid, out_x_r$peakid, nomatch = 0), c("into", "mz", "rt")]),
        as.list(data_r[match(out_x_r$peakid, data_r$tmp_peakid), c("into", "mz", "rt")]))
    }
  }
})


# corPEAKS --------------------------------------------------------------------------------------------------------
test_that("corPEAKS correlates peaks' intensities across samples", {
  ## create dummy table with intensities for 3 peaks across 10 samples:
  ## peaks 1 to 3: same 3 peaks repeated with minor alterations
  peak_1to3_ints <- single_table[single_table$peakgr == 1, "into"]
  samples <- data.frame(peakid = 1:3,
                        into = c(peak_1to3_ints,
                                 peak_1to3_ints - rnorm(3, mean = 100, sd = 20),
                                 peak_1to3_ints - rnorm(3, mean = 100, sd = 20)),
                        run_order = rep(1:3, each = 3))
  samples_cor <- t(utils::combn(unique(samples$peakid), 2, simplify = T))
  weight <- apply(
    samples_cor,
    1,
    FUN = corPEAKS,
    pkg_ints = samples,
    min_samples_n = 0
  )
  expect_true(length(weight) == 3) # 3 peak-peak combinations are posible between 3 peaks
  expect_true(all(weight > 0.8))
  expect_true(class(weight) == "numeric")
  
  ## create dummy table with intensities for 3 peaks across 10 samples:
  ## peaks 1 to 3: same 3 peaks repeated with minor alterations
  ## peak 4: peak has uncorelated intensities
  ## take small into
  samples <- rbind(samples,
                   data.frame(
                     peakid = 4,
                     into = c(200, 250, 300),
                     run_order = 1:3)
                   )
  samples_cor <- t(utils::combn(unique(samples$peakid), 2, simplify = T))
  weight <- apply(
    samples_cor,
    1,
    FUN = corPEAKS,
    pkg_ints = samples,
    min_samples_n = 0
  )
  expect_true(nrow(samples_cor) == 6) # 6 peak-peak combinations are posible between 4 peaks
  expect_true(nrow(samples_cor) == length(weight))
  expect_true(all(weight[which(samples_cor[, 1] == 4 | samples_cor[ , 2] == 4 )] == 0)) # 4th peak is uncorelated
})
