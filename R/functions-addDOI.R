####---- Add DOI peaks to the DB annotations table
## add DATAFILE (doi) to the the TEMPLATE (tmp)
## template could be either the master peak table, or a DATABASE (if building template for the first time)
## returns a template
addDOI <- function(tmp, doi_fname, mz_err, rt_err, bins, add_db = FALSE, db_thrs = 0) {

  ## (A) prepare tmp for grouping with the doi where doi peaks will be added to
  ## get peakgr clusters based on their retention time, order them by complexity
  tmp <- getCLUSTS(dt = tmp)
  tmp <- tmp %>%
    mutate(mz_l = .data$mz - mz_err, mz_h = .data$mz + mz_err, rt_l = .data$rt - rt_err, rt_h = .data$rt + rt_err) %>%
    addCOLS(., c("chemid", "dbid", "dbname", ## if not using db to build a template, add missing columns
                 "doi_peakid", "doi_peakgr", "doi_peakgrcls", "cos")) ## columns for storing intermediate results
  ## dataframe for storing intermediate alignment results
  itmp <- tmp

  ## (B) prepare doi for matching against template
  doi_full <- checkFILE(file = doi_fname)
  doi_full <- getCLUSTS(dt = doi_full)
  ## dataframe for storing intermediate alignment results
  doi <- doi_full %>%
    mutate(mz_l = .data$mz - mz_err, mz_h = .data$mz + mz_err, rt_l = .data$rt - rt_err, rt_h = .data$rt + rt_err) %>%
    mutate(tmp_peakgr = NA, tmp_peakid = NA, chemid = NA, dbid = NA, dbname = NA, cos = NA, added = F) %>%
    select(.data$peakid, .data$mz, .data$mz_l, .data$mz_h, .data$rt, .data$rt_l, .data$rt_h, .data$into, .data$peakgr, .data$peakgrcls, .data$tmp_peakid, .data$tmp_peakgr, .data$chemid, .data$dbid, .data$dbname, .data$cos, .data$added)

  ## initiate progress bar
  pb <- progress_estimated(n = length(unique(doi$peakgrcls)))

  ## (C) for every doi peak, find matching tmp peaks (if any)
  while (any(doi$added == F)) {

    ## get the first peak among the un-annotated peaks
    p <- doi %>% filter(.data$added == F) %>% slice(1) %>% pull(.data$peakid)
    
    target <- doi %>% filter(.data$peakid == p)
    target <- doi %>% filter(.data$peakgrcls == target$peakgrcls)

    ## matching by mz/rt window
    mat <- target %>%
      group_by(.data$peakid) %>%
      do(matchPEAK(p = .data, dt = tmp)) %>%
      ungroup()

    ## extract top matches for every target peakgr
    mattop <- getTOPmatches(mat = mat, target = target, tmp = tmp, bins = bins, add_db = add_db, db_thrs = db_thrs)

    ## add annotated/unmatched DOI peaks to template
    update <- addPEAKS(mattop = mattop, mat = mat, target = target, itmp = itmp, tmp = tmp, doi = doi)

    if (update$itmp %>% filter(!is.na(doi_peakid)) %>% group_by(doi_peakid) %>% summarise(n = n()) %>% filter(n > 1) %>% nrow() > 0) stop("duplicated doi peakid")
    if (update$itmp %>%  group_by(peakid) %>% summarise(n = length(unique(doi_peakgr))) %>% filter(n > 1) %>% nrow() > 0) stop("same peakid to multiple peakgrs")
    if (any(!doi %>% filter(!is.na(tmp_peakid)) %>% pull(peakid) %in% (update$itmp %>% filter(!is.na(doi_peakid)) %>% pull(doi_peakid)))) stop("overwritten doi peakid")
    if (any(is.na(update$itmp$peakgr))) stop("missing peakgr")
    
    itmp <- update$itmp
    doi <- update$doi
  
    ## update progress bar
    pb$tick()$print()

  }

  ## (D) write grouping output
  if (!all(doi$peakid %in% itmp$doi_peakid)) { stop("not all PNOs were assigned to TMP!") }
  message("\n") ## add empty row to align new messages below the progress bar

  doi_full <- right_join(doi_full,
                         itmp %>%
                           filter(!is.na(.data$doi_peakid)) %>%
                           rename(tmp_peakid = .data$peakid, tmp_peakgr = .data$peakgr, peakid = .data$doi_peakid) %>%
                           rename(new_mz = .data$mz, new_rt = .data$rt) %>%
                           select(.data$peakid, .data$tmp_peakid, .data$tmp_peakgr, .data$new_mz, .data$new_rt, .data$chemid, .data$dbid, .data$dbname, .data$cos), by = c("peakid"))
  write.csv(doi_full, file = gsub(".csv", "_aligned.csv", doi_fname), quote = T, row.names = F) ## quote = T to preserve complex DB entries with "-"

  ## (E) Update template
  tmp <- itmp %>%
    ## remove intermediate output columns not needed for next round of alignment
    select(-c(.data$mz_l, .data$mz_h, .data$rt_l, .data$rt_h, .data$peakgrcls, .data$doi_peakid, .data$doi_peakgr, .data$doi_peakgrcls, .data$cos)) %>% 
    ## for any tmp peaks that were matched by multiple doi peakids, retain single, most intense doi peakid
    group_by(.data$peakid) %>% 
    filter(.data$into == max(.data$into)) %>% 
    ungroup()
  
  return(list("tmp" = tmp, "doi" = doi_full))

}
