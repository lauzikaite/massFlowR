## requires columns 'pg', 'pid', rt/mz values
getCLUSTS <- function(dt){
  
  pgc_ids <- data.frame(pg = unique(dt$pg), pgc = NA, stringsAsFactors = F)
  
  ## for every peak group
  for (pgroup in unique(dt$pg)) {
    
    if (is.na(pgc_ids[which(pgc_ids$pg == pgroup),"pgc"])) {
      
      ## assign peak-group-cluster ID
      pgc_id <- ifelse(all(is.na(pgc_ids$pgc)), 1, max(na.omit(pgc_ids$pgc)) + 1)
      
      ## find all other peaks in the same RT region as the peak-group
      peaks <- dt %>%
        filter(pg == pgroup)
      cluster <- dt %>%
        filter(pg %in% (pgc_ids %>% filter(is.na(pgc)) %>% select(pg) %>% pull) ) %>% # only search between peak groups that were not clustered yet
        filter(between(rt, min(peaks$rt), max(peaks$rt))) # use RT region defined by the min and max values
      
      ## add other peaks in the extracted peak-group
      cluster <- bind_rows(cluster,
                       dt %>%
                         filter(pg %in% (cluster %>% distinct(pg) %>% pull(pg))) %>%
                         filter(!pid  %in% (cluster %>% distinct(pid) %>% pull(pid))))
      
      ## update cluster ID assignment table
      pgc_ids[which(pgc_ids$pg %in% unique(cluster$pg)),"pgc"] <- rep(pgc_id, length(unique(cluster$pg)))
      
    } else { next }
    
  }
  
  dt <- full_join(dt, pgc_ids, by = c("pg"))
  
  ## order peak-groups by their pg complexity
  ## and peaks intensity
  dt_c <- dt %>%
    group_by(pg) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    ungroup()
  
  dt <- dt %>%
    arrange(factor(pg, levels = dt_c$pg), pg, pid)
  
  return(dt)
  
}






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
annotateDOI <- function(db, doi, mz_err = 0.01, rt_err = 2, bins = 0.1) {
  
  ## find matches for every DB compound in the 1st study sample
  ## prepare DB annotations for matching
  anno <- db %>% 
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>% 
    arrange(dbid, ) %>% 
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
