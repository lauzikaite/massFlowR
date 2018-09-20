####---- Add DOI peaks to the DB annotations table
## add DATAFILE (doi) to the the TEMPLATE (tmp)
## template could be either the master peak table, or a DATABASE (if building template for the first time)
## returns a template
addDOI <- function(tmp, doi, mz_err, rt_err, bins) {

  ## (A) prepare tmp for grouping with the doi where doi peaks will be added to
  tmp <- getCLUSTS(dt = tmp)
  tmp <- tmp %>%
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err)
  ## build dataframe for storing intermediate alignment results
  itmp <- tmp %>%
    mutate(doi_peakid = NA, doi_peakgr = NA, doi_peakgrcls = NA, cos = NA)

  ## (B) prepare doi for matching against template/db
  doi_full <- getDOI(file = doi)
  ## build dataframe for storing intermediate alignment results
  doi <- doi_full %>%
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
    mutate(tmp_peakgr = NA, tmp_peakid = NA, chemid = NA, dbid = NA, dbname = NA, cos = NA, added = F) %>%
    select(peakid, mz, mz_l, mz_h, rt, rt_l, rt_h, into, peakgr, peakgrcls, tmp_peakid, tmp_peakgr, chemid, dbid, dbname, cos, added)

  ## (C) for every doi peak, find matching tmp peaks (if any)
  while (any(doi$added == F)) {

    ## get the first peak among the un-annotated peaks
    p <- doi %>% filter(added == F) %>% slice(1) %>% pull(peakid)
    print(p)
    target <- doi %>% filter(peakid == p)
    target <- doi %>% filter(peakgrcls == target$peakgrcls)

    ## matching by mz/rt window
    mat <- target %>%
      group_by(peakid) %>%
      do(matchPEAK(p = ., dt = tmp)) %>%
      ungroup()

    ## extract top matches for every target peakgr
    mattop <- getTOPmatches(mat = mat, target = target, tmp = tmp, bins = bins)

    ## add annotated/unmatched DOI peaks to template
    update <- addPEAKS(mattop = mattop, mat = mat, target = target, itmp = itmp, tmp = tmp, doi = doi)

    itmp <- update$itmp
    doi <- update$doi

  }



}
