matchPEAK <- function(p, dt, use_db = T) {

  if (use_db == T) {
    ids <- c("pid", "dbid", "mz", "rt", "into")
  } else {
    ids <- c("cls", pno, pid, mz, rt, into, cid, clid)
  }
  
  mat <- dt %>%
    ## Using regions of both peaks. Regions defined by user mz and rt error ranges
    filter(between(mz_l, p$mz_l, p$mz_h) | between(mz_h, p$mz_l, p$mz_h)) %>%
    filter(between(rt_l, p$rt_l, p$rt_h) | between(rt_h, p$rt_l, p$rt_h)) %>%
    mutate(cls = ifelse(use_db == F, p$cls, NA),
           pno = ifelse(use_db == F, p$pno, NA)) %>% 
    select_(.dots = ids) 
  
  return(mat)

}

####---- Add DOI peaks to the DB annotations table
addDOI <- function(anno, target, mat) {
  
  
  
}




annotateDOI <- function(db, doi, mz_err = 0.01, rt_err = 2, bins = 0.1) {
  
  ## find matches for every DB compound in the 1st study sample
  ## 'pid'numbering starts in DB anno
  ## prepare DB annotations gor grouping 
  anno <- db %>% 
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>% 
    arrange(dbid) %>% 
    mutate(pid = row_number()) %>% #, comp = NA, cls = NA, cid = NA, clid = NA, cos = NA) %>% 
    select(-pno)

  ## prepare 1st study sample for grouping
  doi_full <- read.csv(doi, stringsAsFactors = F)
  doi <- doi_full %>%
    select(pno, mz, rt, into, comp) %>%
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err,
           cls = NA)
  
  doi_peaks <- doi %>%
    mutate(dbid = NA, matched = F) %>% 
    select(pno, dbid, comp, matched)
  
  
  while (any(doi_peaks$matched == F)) {
  
    p <- doi_peaks %>% filter(matched == F) %>% slice(1) %>% pull(pno)
    target <- doi %>% filter(pno == p)
    target <- doi %>% filter(comp == target$comp)
    
    ## matching by mz/er window
    mat <- target %>%
      group_by(pno) %>%
      do(matchPEAK(p = ., dt = anno, use_db = T)) %>% 
      ungroup()
    
    ## add DOI peaks to DB annotations table
    out <- addDOI()
    anno <- out$anno
    doi_peaks <- out$doi_peaks
      
    }
    
   
  
}
