context("addDOI")

test_that("Selection of top matches via compareCLUSTERS() is correct", {

  ## sample table consists of 3 target peakgrs, each matched by the same 3 template peakgrs
  cos_list <- c(seq(from = 0.1, to = 1, length.out = 8), 1)
  matcos <- data.frame(target_peakgr = c(rep(1,3), rep(2,3), rep(3,3)),
                       target_peakgrcls = c(rep(1, 9)),
                       peakgr = c(rep(c(1,2,3), 3)),
                       peakgrcls = c(rep(1, 9)),
                       chemid = c(rep(NA, 9)),
                       cos = cos_list,
                       stringsAsFactors = F)
  expect_error(top <- compareCLUSTERS(matcos = matcos, add_db = F), "identical cosines were found!")

  ## 1-1 is top pair by rank
  ## 2-2 is second best for target peakgr (choosing from what is left unassigned)
  ## 3-3 is third best for target peakgr (choosing from what is left unassigned)
  cos_list <- seq(from = 0.1, to = 1, length.out = 9)
  cos_list <- cos_list[c(9,6,3, 8,5,2, 7,4,1)]
  matcos <- matcos %>%
    mutate(cos = cos_list)

  top <- compareCLUSTERS(matcos = matcos, add_db = F)
  expected_top <- data.frame(target_peakgr = c(1,2,3),
                             target_peakgrcls = rep(1,3),
                             peakgr = c(1,2,3),
                             peakgrcls = rep(1,3),
                             chemid = rep(NA,3),
                             cos = c(1,0.55,0.1),
                             top = rep(TRUE,3))
  expect_equal(top, expected_top)

  })
