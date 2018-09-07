## requires columns 'peakg', 'peakid', rt/mz values
getCLUSTS <- function(dt){

  ## data frame to store assigned peak-group-clusters' ids
  cls <- data.frame(peakgr= unique(dt$peakgr), peakgrcls = NA, stringsAsFactors = F)

  ## for every peak group
  for (pg in unique(dt$peakgr)) {

    ## if peak-group was not assigned to a cluster yet
    if (is.na(cls[which(cls$peakgr == pg),"peakgrcls"])) {

      ## get unique peak-group-cluster id
      id <- ifelse(all(is.na(cls$peakgrcls)), 1, max(na.omit(cls$peakgrcls)) + 1)

      ## find all other peaks in the same RT region as the peak-group
      peaks <- dt %>%
        filter(peakgr == pg)
      cluster <- dt %>%
        filter(peakgr %in% (cls %>% filter(is.na(peakgrcls)) %>% select(peakgr) %>% pull) ) %>% # only search between peak groups that were not clustered yet
        filter(between(rt, min(peaks$rt), max(peaks$rt))) # use RT region defined by the min and max values

      ## add other peaks in the extracted peak-group
      cluster <- bind_rows(cluster,
                       dt %>%
                         filter(peakgr %in% (cluster %>% distinct(peakgr) %>% pull(peakgr))) %>%
                         filter(!peakid  %in% (cluster %>% distinct(peakid) %>% pull(peakid))))

      ## update cluster ID assignment table
      cls[which(cls$peakgr %in% unique(cluster$peakgr)),"peakgrcls"] <- rep(id, length(unique(cluster$peakgr)))

    } else { next }

  }

  dt <- full_join(dt, cls, by = c("peakgr"))

  ## order peak-groups by their peakgrcomplexity
  ## and peaks intensity
  dt_c <- dt %>%
    group_by(peakgr) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    ungroup()

  dt <- dt %>%
    arrange(factor(peakgr, levels = dt_c$peakgr), peakgr, peakid)

  return(dt)

}





## find matches for a peak from one dataset in another
## use peak's and dataset's mz and rt regions only
matchPEAK <- function(p, dt) {
  mat <- dt %>%
    ## Using regions of both peaks. Regions defined by user mz and rt error ranges
    filter(between(mz_l, p$mz_l, p$mz_h) | between(mz_h, p$mz_l, p$mz_h)) %>%
    filter(between(rt_l, p$rt_l, p$rt_h) | between(rt_h, p$rt_l, p$rt_h)) %>%
    mutate(target_peakid = p$peakid, target_peakgr = p$peakgr) %>%
    select(target_peakid, target_peakgr, names(dt))

  return(mat)

}

####---- Add DOI peaks to the DB annotations table
annotateDOI <- function(db, doi, mz_err = 0.01, rt_err = 2, bins = 0.1) {

  ## find matches for every DB compound in the 1st study sample
  ## prepare DB annotations for matching
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
