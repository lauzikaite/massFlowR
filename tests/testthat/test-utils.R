context("utils functions")

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
  expect_true(all_equal(test_pks_rd %>% select(-peakid),
                        pks_clean))
  
  ## removes duplicated peakids (as in addDOI)
  # create peak table with duplicated peakids
  test_pks_dup <- test_pks_rd %>%
    bind_rows(., test_pks_rd %>% filter(peakid %in% c(1, 10, 20)))
  # cleanup peak table
  pks_unique <- unique(test_pks_dup[, c("mz", "rt")])
  pks_clean <- lapply(
    1:nrow(pks_unique),
    FUN = cleanPEAKS,
    dt_unique = pks_unique,
    dt = test_pks_dup
  )
  pks_clean <- do.call("rbindCLEAN", pks_clean)
  expect_true(any(duplicated(test_pks_dup$peakid)))
  expect_true(all_equal(test_pks_rd, pks_clean))
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
  random_peak <- sample(test_pks_rd$peakid, size = 1)
  random_peak_scpos <-
    c(test_pks_rd[random_peak, "scpos"] - 1, test_pks_rd[random_peak, "scpos"] + 1)
  random_peak_co <- test_pks_rd[
    which(dplyr::between(test_pks_rd$scpos,  random_peak_scpos[1],  random_peak_scpos[2])),
    "peakid"]
  
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
  expect_true(all(apply(random_peak_pairs, 1, function(pair) {
    nrow(
      getCORmat_out %>%
        filter(from == pair[1] | to == pair[1]) %>%
        filter(from == pair[[2]] | to == pair[[2]])
    )
  }) == 1))

  })

# buildGRAPH ------------------------------------------------------------------------------------------------------
test_that("buildGRAPH", {
  
  ## use indeces of co-eluting peaks for a specific peak no 21
  peak <- 21
  peak_scpos <-
    c(test_pks_rd[peak, "scpos"] - 1, test_pks_rd[peak, "scpos"] + 1)
  peak_co <-
    test_pks_rd[which(dplyr::between(test_pks_rd$scpos,  peak_scpos[1],  peak_scpos[2])),
                "peakid"]
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