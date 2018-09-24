####---- Add DOI peaks to the DB annotations table
## add DATAFILE (doi) to the the TEMPLATE (tmp)
## template could be either the master peak table, or a DATABASE (if building template for the first time)
## returns a template
addDOI <- function(tmp, doi_fname, mz_err, rt_err, bins, add_db = FALSE, db_thrs = 0) {

  ## (A) prepare tmp for grouping with the doi where doi peaks will be added to
  tmp <- getCLUSTS(dt = tmp)
  tmp <- tmp %>%
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
    addCOLS(., c("chemid", "dbid", "dbname")) %>%  ## if not using db to build a template, add missing columns
    mutate(doi_peakid = NA, doi_peakgr = NA, doi_peakgrcls = NA, cos = NA) ## if tmp has been aligned with doi before, remove old metadata

  ## build dataframe for storing intermediate alignment results
  itmp <- tmp

  ## (B) prepare doi for matching against template/db
  doi_full <- getDOI(file = doi_fname)
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
    mattop <- getTOPmatches(mat = mat, target = target, tmp = tmp, bins = bins, add_db = add_db, db_thrs = db_thrs)

    ## add annotated/unmatched DOI peaks to template
    update <- addPEAKS(mattop = mattop, mat = mat, target = target, itmp = itmp, tmp = tmp, doi = doi)

    itmp <- update$itmp
    doi <- update$doi

  }

  ## (D) write grouping output
  if (!all(doi$peakid %in% itmp$doi_peakid)) { stop("not all PNOs were assigned to TMP!") }

  doi_full <- right_join(doi_full,
                         itmp %>%
                           filter(!is.na(doi_peakid)) %>%
                           rename(tmp_peakid = peakid, tmp_peakgr = peakgr, peakid = doi_peakid) %>%
                           rename(new_mz = mz, new_rt = rt) %>%
                           select(peakid, tmp_peakid, tmp_peakgr, new_mz, new_rt, chemid, dbid, dbname, cos), by = c("peakid"))
  write.csv(doi_full, file = gsub(".csv", "_aligned.csv", doi_fname), quote = T, row.names = F) ## quote = T to preserve complex DB entries with "-"

  ## Update template
  tmp <- itmp %>%
    select(-peakgrcls) ## must be removed, otherwise next addDOI call with have two peakgrcls columns

  return(list("tmp" = tmp, "doi" = doi_full))

}
