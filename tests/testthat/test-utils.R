context("utils functions")

# validFILE --------------------------------------------------------------------------------------------------------
test_that("validFILE", {
  # original CDF
  expect_true(validFILE(test_fname)) 
  # incorrect CDF filename
  expect_equal(validFILE(gsub(".CDF", ".xxx", test_fname)),
               paste0("File doesn't exist: ", gsub(".CDF", ".xxx", test_fname)))
  ## csv is a wrong filetype
  expect_equal(validFILE(meta_fname), 
               paste0("File format is unsupported: ", meta_fname))
})

# readDATA --------------------------------------------------------------------------------------------------------
test_that("readDATA", {
  raw <- readDATA(test_fname)
  expect_true(class(raw) ==  "OnDiskMSnExp")
})

# cleanPEAKS ------------------------------------------------------------------------------------------------------
test_that("cleanPEAKS removes duplicated peaks correctly", {
  ## removes duplicates in mz/rt (as in pickPEAKS)
  # use un-cleaned, original centwave output
  pks_unique <- unique(test_pks[, c("mz", "rt")])
  pks_clean <- lapply(
    1:nrow(pks_unique),
    FUN = cleanPEAKS,
    dt_unique = pks_unique,
    dt = test_pks
  )
  pks_clean <- do.call("rbindCLEAN", pks_clean)
  pks_clean <- pks_clean[order(pks_clean$into, decreasing = T), ]
  pks_clean$peakid <- 1:nrow(pks_clean)
  
  expect_true(nrow(pks_clean) == nrow(test_pks_rd))
  expect_true(all(sapply(1:nrow(pks_clean), function(x) {
    pks_clean$mz[x] == test_pks_rd$mz[x]
  })))
  expect_true(all(sapply(1:nrow(pks_clean), function(x) {
    pks_clean$rt[x] == test_pks_rd$rt[x]
  })))
  
  ## removes duplicated peakids (as in addDOI)
  # create peak table with duplicated peaks 1, 10 and 20
  test_pks_dup <- rbind(test_pks_rd,
                        test_pks_rd[which(test_pks_rd$peakid %in% c(1, 10, 20)), ])
  pks_unique <- unique(test_pks_dup[, c("mz", "rt")])
  pks_clean <- lapply(
    1:nrow(pks_unique),
    FUN = cleanPEAKS,
    dt_unique = pks_unique,
    dt = test_pks_dup
  )
  pks_clean <- do.call("rbindCLEAN", pks_clean)
  expect_true(nrow(pks_clean) == nrow(test_pks_rd))
  expect_true(any(duplicated(test_pks_dup$peakid)))
  expect_true(all(sapply(1:nrow(pks_clean), function(x) {
    pks_clean$mz[x] == test_pks_rd$mz[x]
  })))
  expect_true(all(sapply(1:nrow(pks_clean), function(x) {
    pks_clean$rt[x] == test_pks_rd$rt[x]
  })))
})

# scaleEDGES ------------------------------------------------------------------------------------------------------
test_that("scaleEDGES with default values", {
  cor_cofs <- c(
    0.97066239157605,
    0.999263924873744,
    0.968474483893108,
    0.980828448489916,
    0.991226243530123,
    0.972021809991886,
    0.998102273445145,
    0.946544580264379,
    0.950924151963622,
    0.969105442905998
  )
  ## using default values
  scaleEDGES_out <- scaleEDGES(x = cor_cofs)
  expect_equal(which(cor_cofs == min(cor_cofs)), which(scaleEDGES_out == min(scaleEDGES_out)))
  expect_equal(min(scaleEDGES_out), 0.01)
  expect_equal(which(cor_cofs == max(cor_cofs)), which(scaleEDGES_out == max(scaleEDGES_out)))
  expect_equal(max(scaleEDGES_out), 10)
})

# rbindCLEAN ------------------------------------------------------------------------------------------------------
test_that("rbindCLEAN", {
  ## make a list of tables
  # replicate same table ntimes
  ntimes <- 10
  list_of_pks <- lapply(1:ntimes, function(dt) {
    test_pks_rd
  })
  rbindCLEAN_out <- do.call("rbindCLEAN", list_of_pks)
  expect_true(all(rownames(rbindCLEAN_out) == 1:nrow(rbindCLEAN_out)))
  expect_equal(length(unique(rbindCLEAN_out$peakid)), max(test_pks_rd$peakid))
  expect_equal(length(which(rbindCLEAN_out$peakid == 1)), ntimes)
})

# getCORmat -------------------------------------------------------------------------------------------------------
test_that("getCORmat", {
  ## use indeces of co-eluting peaks for a random peak
  # random_peak <- sample(test_pks_rd$peakid, size = 1)
  random_peak <- 91
  random_peak_scpos <-
    c(test_pks_rd[random_peak, "scpos"] - 1, test_pks_rd[random_peak, "scpos"] + 1)
  random_peak_co <-
    test_pks_rd[which(test_pks_rd$scpos >= random_peak_scpos[1] &
                        test_pks_rd$scpos <= random_peak_scpos[2]), "peakid"]
  
  ## basic matching
  getCORmat_out <- massFlowR:::getCORmat(ind = random_peak_co)
  expect_true(all(names(getCORmat_out) == c("from", "to")))
  expect_true(all(random_peak_co %in% unique(c(getCORmat_out$from, getCORmat_out$to))))
  
  ## get unique co-eluting peak-peak combinations
  random_peak_pairs <- lapply(random_peak_co, function(peak) {
    other <- random_peak_co[which(random_peak_co != peak)]
    pairs <- lapply(other, function(o) {
      c(peak, o)
    })
    do.call("rbind", pairs)
  })
  random_peak_pairs <- do.call("rbind", random_peak_pairs)
  random_peak_pairs_out <- lapply(1:nrow(random_peak_pairs), function(x) {
    pair <- random_peak_pairs[x, ]
    pair_mat <- getCORmat_out[which((getCORmat_out$from == pair[1] |
                                       getCORmat_out$to == pair[1]) &
                                      (getCORmat_out$from == pair[2] |
                                         getCORmat_out$to == pair[2])), ]
    nrow(pair_mat)
  })
  expect_true(all(random_peak_pairs_out == 1))
  })

# buildGRAPH ------------------------------------------------------------------------------------------------------
test_that("buildGRAPH", {
  
  ## use indeces of co-eluting peaks for a specific peak no 21
  peak <- 21
  peak_scpos <-
    c(test_pks_rd[peak, "scpos"] - 1, test_pks_rd[peak, "scpos"] + 1)
  peak_co <- test_pks_rd[which(test_pks_rd$scpos >= peak_scpos[1] &
                                 test_pks_rd$scpos <= peak_scpos[2]),  "peakid"]
  peak_cormat <- getCORmat(ind = peak_co)
  peak_cormat$weight <-
    apply(peak_cormat, 1, FUN = corEIC, eic = test_eic_rd)
  buildGRAPH_out <- buildGRAPH(pkg_cor = peak_cormat,
                               cor_thr = 0.95,
                               plot = FALSE)
  expect_true(class(buildGRAPH_out) == "membership")
  expect_true(all(peak_co %in% names(buildGRAPH_out)))
  expect_true(all(buildGRAPH_out >= 1))

})
