addGROUPED <- function(component, tmat, tmp, mz_err, rt_err, doi_peaks) {

  ## take the TMP CID, which was assigned to the component
  component <- component %>%
    filter(!is.na(topCOMP))

  ####---- (1) extract matches-target table for the COMP-CID combination
  componentmat <- tmat %>%
    filter(
      (target == T & comp == component$comp)|
        (target == F & comp == component$comp & cid == component$cid)|
        (target == F & is.na(comp) & cid == component$cid)) %>%
    mutate(
      ## add matching evaluation information
      cid = component$cid,
      cos = unlist(component$cos),
      tmpCLS_order = component$tmpCLS_order,
      targetCLS_sc = component$targetCLS_sc,
      tmpCLS_sc = component$tmpCLS_sc)

  ####---- (2) merge TMP peaks that match to the same KEY

  ## (2A) list all KEYS that the tmp peak is matched to
  cmat <- componentmat %>%
    group_by(target, pid) %>% # pid is either from TMP or DOI
    arrange(pid) %>% # order by pid (intensity of the feature)
    mutate(key_paste = ifelse(is.na(key), NA, paste0(key, collapse = ",")))

  ## (2B) for each KEY that matches multiple tmp peaks, find the closest peak
  cmat <- cmat %>%
    group_by(target, key) %>%
    do(selectPEAK(t = ., cmat = cmat)) %>%
    ungroup()

  ## (2C) for each KEY, average the mz/rt across all matches
  cmat <- cmat %>%
    group_by(key) %>%
    mutate(
      mz_final = ifelse(is.na(key), mz,
                        ifelse(any(target == FALSE), median(mz), mz)),
      rt_final = ifelse(is.na(key), rt,
                        ifelse(any(target == FALSE), median(rt), rt))) %>%
    ungroup() %>%
    mutate(
      cid,
      mz = mz_final, rt = rt_final, # change mz/rt to the median of mz/rt across all matches for that peak
      ## add matched DOI component info
      cls = tmat %>% filter(target == T) %>% distinct(clid) %>% pull())

  ## (2D) for each PEAK, if more than one KEY matched, average mz/rt across all matched KEYS and retain a single row
  cmat <- cmat %>%
    group_by(target, pid) %>%
    mutate(
      mz_final = if(all(target == F)) { ifelse(is.na(key), mz, median(mz)) } else { mz },
      rt_final = if(all(target == F)) { ifelse(is.na(key), rt, median(rt)) } else { rt }) %>%
    filter((target == FALSE & row_number() == 1) | target == TRUE) %>%
    ungroup() %>%
    mutate(
      cid,
      ## change mz/rt to the median of mz/rt across all matches for that peak
      mz = mz_final, rt = rt_final,
      mz_l = mz_final - mz_err, mz_h = mz_final + mz_err, rt_l = rt_final - rt_err, rt_h = rt_final + rt_err)

  ####---- (3)  update matched TMP peaks in the TMP
  utmp <- cmat %>%
    ## since each key-match combination now has a single median mz/rt value, a single feature can be retained for the key
    filter(target == FALSE) %>%
    mutate(tmp = TRUE) %>%
    mutate(clid = component$clid) %>%
    select(pid, mz, rt, into, cid, clid, mz_l, mz_h, rt_l, rt_h, tmp, pno, key = key_paste, key_max, comp, cls, cos, tmpCLS_order, targetCLS_sc, tmpCLS_sc)

  ####---- (A4) add unmatched DOI target peaks to TMP as new peaks under the same CID

  ## retain only those target peaks that did not have a match
  un_cmat <- componentmat %>%
    group_by(key) %>%
    filter(ifelse(all(target == T), T, F)) %>%
    ungroup()

  if(nrow(un_cmat) > 0) {

    ## generate new PIDs for these DOI peaks
    pid_id <- (max(tmp$pid) + 1):(max(tmp$pid) + nrow(un_cmat))

    un_cmat <- un_cmat %>%
      mutate(
        tmp = FALSE,
        pno = pid,
        pid = pid_id,
        cls = clid, # DOI cluster id
        clid = component$clid) %>%
      select(pid, mz, rt, into, cid, clid, mz_l, mz_h, rt_l, rt_h, tmp, pno, key, key_max, pno, comp, cls, cos, tmpCLS_order, targetCLS_sc, tmpCLS_sc)

    ## add all matches into one df
    utmp <- bind_rows(utmp, un_cmat)
  }

  ## assign CIDs for the checked DOI peaks
  doi_peaks[which(doi_peaks$pno %in% (cmat %>% filter(target == T) %>% pull(pno))),] <- utmp %>%
    filter(pno %in% (cmat %>% filter(target == T) %>% pull(pno))) %>%
    select(pno, pid, cid, clid)

  return(list("doi_peaks" = doi_peaks, "utmp" = utmp))

}

addUNGROUPED <- function(component, tmat, tmp, mz_err, rt_err, doi_peaks) {

  ## since DOI component was not asssigned, it doesn't matter how many CIDs were tested for, and thus can retain just one row
  component <- component %>%
    slice(1)

  ## extract matches-target table for the single target's component
  cmat <- tmat %>%
    filter(target == T & comp == unique(component$comp))

  ## assign new CID and PID for the unmatched target component
  cid_id <- max(tmp$cid) + 1
  pid_id <- (max(tmp$pid) + 1):(max(tmp$pid) + nrow(cmat))

  utmp <- cmat %>%
    mutate(pid = pid_id,
           clid = NA, # temporalily assign cluster ID with NA, since ids will be generated again at the end of tmp grouping round
           mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err,
           tmp = FALSE,
           ## since keys are not pasted, numeric format would conflict with other scenarios
           key = as.character(key), key_max = as.character(key_max),
           ## DOI's info
           comp = cid,
           cos = NA,
           ## TMP's info comes last since this is erasing
           cid = cid_id,
           tmpCLS_order = NA,
           targetCLS_sc = component$targetCLS_sc,
           tmpCLS_sc = NA) %>%
    select(pid, mz, rt, into, cid, clid, mz_l, mz_h, rt_l, rt_h, tmp, pno, key, key_max, comp, cls, cos, tmpCLS_order, targetCLS_sc, tmpCLS_sc)

  # utmp <- addUNGROUPEDnew(cmat = cmat, tmp = tmp, component = component, mz_err = mz_err, rt_err = rt_err)

  if( nrow(cmat) != nrow(utmp)) { stop(paste("un-assigned DOI component peaks are not added to TMP correctly! Check p:", p)) }

  ## assign CIDs for the checked DOI peaks
  doi_peaks[which(doi_peaks$pno %in% (cmat %>% filter(target == T) %>% pull(pno))),] <- utmp %>%
    filter(pno %in% (cmat %>% filter(target == T) %>% pull(pno))) %>%
    select(pno, pid, cid, clid)

  return(list("doi_peaks" = doi_peaks, "utmp" = utmp))

}


