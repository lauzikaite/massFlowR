
####---- standartise the template/doi before each new round of grouping
getCLUSTS <- function(dt){

  ## build cluster id matching table
  cls_ids <- data.frame(cid = unique(dt$cid), clid = NA)

  for (cid_id in unique(unique(dt$cid))) {

    ## if component was not grouped yet into a cluster
    if (is.na(cls_ids[which(cls_ids$cid == cid_id),"clid"])) {

      ## assign cluster ID
      cls_id <- ifelse(all(is.na(cls_ids$clid)), 1, max(na.omit(cls_ids$clid)) + 1)

      ## extract component's peaks
      cmp <- dt %>%
        filter(cid == cid_id)

      ## find all other peaks in the same RT region
      mat <- dt %>%
        filter(cid %in% (cls_ids %>% filter(is.na(clid)) %>% select(cid) %>% pull) ) %>% # only search between components that were not clustered yet
        filter(between(rt, min(cmp$rt), max(cmp$rt))) # use RT region defined by the min and max values for all component's peaks

      ## add other peaks in the extracted components
      mat <- bind_rows(mat,
                       dt %>%
                         filter(cid %in% (mat %>% distinct(cid) %>% pull(cid))) %>%
                         filter(!pid  %in% (mat %>% distinct(pid) %>% pull(pid))))

      ## update cluster ID assignment table
      cls_ids[which(cls_ids$cid %in% unique(mat$cid)),"clid"] <- rep(cls_id, length(unique(mat$cid)))

    } else { next }

  }

  dt <- full_join(dt, cls_ids, by = c("cid"))

  ## order peaks by their components' complexity (components with more peaks go first),
  ## and peaks intensity (more intense peaks in the component go first)
  dt_c <- dt %>%
    group_by(cid) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    ungroup()

  dt <- dt %>%
    arrange(factor(cid, levels = dt_c$cid), cid, pid)

  return(dt)

}


####---- getMATCH(): makes long df wit single-row per feature, either target, or match. Column 'target' indicates where the peak comes from
getMATCH <- function(t, tmp) {

  mat <- tmp %>%
    ## Using regions of both peaks. Regions defined by user mz and rt error ranges
    filter(between(mz_l, t$mz_l, t$mz_h) | between(mz_h, t$mz_l, t$mz_h)) %>%
    filter(between(rt_l, t$rt_l, t$rt_h) | between(rt_h, t$rt_l, t$rt_h)) %>%

    ## only match to peaks that were already in the tmp and did not come from DOI
    ## this is required as template is iterated over the same DOI peak table
    filter(is.na(tmp) | tmp == T) %>%
    mutate(comp = t$comp, cls = t$cls, pno = t$pno, key = t$key, key_max = t$key_max, target = FALSE) %>%
    select(comp, cls, key, key_max, pno, mz, rt, into, pid, cid, clid, target)

  if (nrow(mat) == 0) {

    mat <- t %>%
      mutate(pid = pno, target = TRUE, cid = comp,  clid = cls) %>%
      select(comp, cls, key, key_max, pno, mz, rt, into, pid, cid, clid, target)


  } else {

    mat <- bind_rows(
      # assign original DOI component to mat cid
      t %>%
        mutate(pid = pno, target = TRUE, cid = comp,  clid = cls) %>%
        select(comp, cls, key, key_max, pno, mz, rt, into, pid, cid, clid, target),
      mat)
  }
  return(mat)

}

####---- scale vector to unit length
scaleVEC <- function(x) {
  x / ( sqrt(sum(x * x)) )
}

####---- cosine estimation
getCOS <- function(breaks, target_vec, vector, cluster){

  vector_vec <- full_join(breaks,
                          vector,
                          by = c("bin")) %>%
    mutate(into = ifelse(is.na(into), 0, into)) %>%
    select(pid, breaks, bin, into)

  ## remove duplicated vector PIDs, that are vectored to multiple key - RETAIN ORIGINAL ROW ORDER BY BIN NUMBER
  vector_vec <- bind_rows(
    ## retain all target_vec entries
    vector_vec %>% filter(is.na(pid)),
    ## remove duplicating PIDs and/or duplicating bins from different PIDS
    vector_vec %>% filter(!is.na(pid)) %>% group_by(bin) %>% filter(row_number(pid) == 1)
    # vector_vec %>% filter(!is.na(pid)) %>% group_by(pid) %>% filter(row_number(pid) == 1)
  ) %>% arrange(bin)

  ## scale intensities
  vector_vec <- scaleVEC(vector_vec$into)

  ## find the cosine of the angle
  cos_angle <- (sum(target_vec * vector_vec) )  / ( (sqrt(sum(target_vec * target_vec)))  * ( sqrt(sum(vector_vec * vector_vec)) ) )
  return(cos_angle)
}

####---- get matching scenario
getSCEN <- function(tmat) {

  ## (1) check for CID multiplicicity - do same CID get matched by multiple COMP?
  multi <- tmat %>%
    filter(target == F) %>%
    group_by(clid, cid) %>%
    summarise(target_comp = paste0(unique(comp), collapse = ","),
              target_comp_n = length(unique(comp)))

  if(any(multi$target_comp_n > 1)) {

    ## Scenario C: clusters are compared together
    "C"

  } else {

    ## (2) check for COMP multiplicity - do same COMP match multiple CID?
    multi <- tmat %>%
      filter(target == F) %>%
      group_by(comp) %>%
      summarise(match_cid = paste0(unique(cid), collapse = ","),
                match_cid_n = length(unique(cid)))

    if( any(multi$match_cid_n > 1)) {

      ## Scenario B: decision on CID/CLS
      "B"

    } else {

      ## Scenario A: one-to-one matching, no decision needed
      "A"
    }
  }
}

####---- Compare target COMP and all matched TMP CIDs
compareCOMPS <- function(t, tmat, bins) {

  ## if this component have any matches
  cmat <- tmat %>% filter(target == F, comp %in% t$comp) ## matches by mz/rt

  if(nrow(cmat) > 0) {

    cmat <-  bind_rows(t, # target peaks
                       cmat, # tmp peaks matched by target's mz/rt
                       tmat %>% filter(target == F, is.na(comp), cid %in% cmat$cid) # tmp peaks matched by cid only
    )

    breaks <- data.frame(breaks = seq(from = min(cmat$mz), to = max(cmat$mz) , by = bins)) %>%
      mutate(bin = row_number())
    cmat <- cmat %>%
      mutate(bin = findInterval(mz, breaks$breaks))

    ## build target vector of the target spectrum
    target_vec <- full_join(breaks,
                            cmat %>%
                              filter(target == TRUE) %>%
                              ## remove duplicating duplicating bins from different PIDS
                              group_by(bin) %>% filter(row_number(pid) == 1) %>%
                              select(bin, into),
                            by = c("bin")) %>%
      mutate(into = ifelse(is.na(into), 0, into))

    ## scale vector to unit length of 1
    target_vec <- scaleVEC(target_vec$into)

    ## get the cosine of the angle between vectors representing spectra
    cos <- cmat %>%
      filter(target == FALSE) %>%
      group_by(clid, cid) %>%
      ## do returns a dataframe with two columns: first is the label, second is a list with outcomes of the function
      do(cos = getCOS(breaks = breaks, target_vec = target_vec, vector = .)) %>%
      ungroup()


  } else {

    cos <- data.frame(clid = NA, cid = NA, cos = NA)
  }

  return(cos)
}

####---- Compare target CLS and all matched TMP CLSs, extract top-matched TMP CIDs
compareCLS <- function(matcos, tmat) {

  ## extract TOP matches for each target COMP (if there are any)
  mattop <- matcos %>%
    group_by(cid) %>%
    mutate(multiCIDn = n()) %>%

    ## (1) for each CID, extract only one COMP, based on cosine
    mutate(topCID =
             ## retain COMPS with no CID
             if(all(!is.na(cid))) {
               ## top from multiple COMP
               ifelse(multiCIDn > 1 & cos == max(unlist(cos)), TRUE,
                      ## top from single COMP
                      ifelse(multiCIDn == 1 & cos == max(unlist(cos)), TRUE, NA)) } else { NA}) %>%

    ## (2) for each COMP, retain only one CID
    group_by(topCID, comp) %>%
    mutate(topCIDn = sum(topCID == TRUE)) %>%
    ## select topCOMP from single/multiple topCID
    mutate(topCOMP = if(all(!is.na(cid))) { ifelse(topCIDn > 0 & topCID == TRUE & cos == max(unlist(cos)), TRUE, NA) } else { NA }) %>%
    group_by(comp) %>%
    ## select topCOMP when no topCID is available
    mutate(topCOMP = ifelse(all(is.na(topCOMP)) & !is.na(cid) & multiCIDn == 1 & n() == 1, TRUE, topCOMP)) %>%
    ungroup()

  ## add those TMP cids, which did not have any match, but are in the same CLS
  mattop <- bind_rows(mattop,
                      tmat %>%
                        filter(target == F, !cid %in% matcos$cid) %>%
                        group_by(cid, clid) %>%
                        summarise(comp = NA, cos = NA))

  ## get order of components in clusters based on RT
  o_target <- tmat %>%
    filter(target == T) %>%
    group_by(comp) %>%
    summarise(rtmed = median(rt)) %>%
    arrange(rtmed) %>%
    ungroup() %>%
    mutate(order = row_number())

  o_tmp <- tmat %>%
    filter(target == F) %>%
    group_by(clid, cid) %>% # column 'clid' contains cluster ids for tmp/doi, depending on row
    summarise(rtmed = median(rt)) %>%
    arrange(rtmed) %>%
    group_by(clid) %>%
    mutate(order = row_number()) %>%
    ungroup()

  ## add columns for components order in the clusters of target and tmp
  mattop <-
    full_join(
      full_join(
        mattop,
        o_target %>% select(comp, target_order = order),
        by = c("comp")),
      o_tmp %>% select(cid, tmp_order = order),
      by = c("cid")
    )

  ## cluster summaries
  matcls <- mattop %>%
    group_by(clid) %>%
    arrange(target_order) %>%
    filter(topCOMP == T) %>%
    summarise(tmpCLS_order = all(tmp_order == cummax(tmp_order)), # TRUE if every matched comp/cid pair stands in the same position in the cluster
              targetCLS_sc = sum(!is.na(topCOMP))/ max(o_target$order), # ratio of matched comps to cid, to total cids in the cluster
              tmpCLS_sc =  sum(!is.na(topCOMP))/ max(tmp_order)) # ratio of matched comps to cid, to total comps in the cluster
  ## add summaries
  mattop <- full_join(mattop %>% select(comp, clid, cid, cos, topCOMP),
                      matcls, by = c("clid"))

  ## for each COMP, retain only the top selected CID (if none, retain one row per COMP)
  mattop <- mattop %>%
    group_by(comp) %>%
    filter(topCOMP == T | ifelse(all(is.na(topCOMP)), is.na(topCOMP), NA)) %>%
    ungroup()

  return(mattop)

}


####---- Perform matching comparison and template update
runSCEN <- function(tmat, scen, tmp, bins, mz_err, rt_err, doi_peaks) {

  ## (1) add additional peaks from the matched TMP clusters
  tmat <- bind_rows(tmat,
                    tmp %>%
                      ## ! only use peaks that were originally in the tmp: filter(is.na(tmp) | tmp == T)
                      filter( (is.na(tmp) | tmp == T)) %>%
                      ## add additional tmp cls peaks
                      filter(clid %in% (tmat %>% filter(target == FALSE, !is.na(clid)) %>% distinct(clid) %>% pull(clid))) %>%
                      filter(!pid  %in% (tmat %>% filter(target == FALSE, !is.na(clid)) %>% distinct(pid) %>% pull(pid))) %>%
                      mutate(target = FALSE))

  ####--- (A|B) simple component-by-component scenario

  if(any(scen == "A" | scen == "B")) {

    ## (1) compare all target COMPS against TMP, one-by-one
    mattop <- tmat %>%
      filter(target == T) %>%
      group_by(comp) %>%
      do(compareCOMPS(t = ., tmat = tmat, bins = bins)) %>%
      ## if more than one CID is matched, take the one with highest cosine, NA is used instead of FALSE since in compareCLS() NA is generated
      mutate(topCOMP = ifelse(is.na(cid), NA, T)) %>%
      filter(if (all(!is.na(topCOMP))) { cos == max(unlist(cos)) } else { is.na(topCOMP) }) %>%
      ungroup() %>%
      ## add this for updateTMP to work correctly
      mutate(tmpCLS_order = NA, targetCLS_sc = NA, tmpCLS_sc = NA)


    ####---- (C) complex full-cluster matching scenario
  } else {

    ## (1) check similarity of all target COMPS against TMP, one-by-one
    matcos <- tmat %>%
      filter(target == T) %>%
      group_by(comp) %>%
      do(compareCOMPS(t = ., tmat = tmat, bins = bins))

    ## (2) compare CLS, extract top matches and estimate match scores
    mattop <- compareCLS(matcos = matcos, tmat = tmat)

    if (!all(tmat %>% filter(target == T) %>% distinct(comp) %>% pull() %in% mattop$comp)) { stop(paste0("compareCLS does not add all DOI components! Check p:", p)) }

  }

  ####---- universal part for updating TMP
  # update <- updateTMP(mattop = mattop, tmat = tmat, doi_peaks = doi_peaks, tmp = tmp, mz_err = mz_err, rt_err = rt_err)
  update <- addCOMPS(mattop = mattop, tmat = tmat, doi_peaks = doi_peaks, tmp = tmp, mz_err = mz_err, rt_err = rt_err)

  return(list("doi_peaks" = update$doi_peaks, "tmp" = update$tmp))

}

####---- Update template with selected matches
# updateTMP <- function(mattop, tmat, doi_peaks, tmp, mz_err, rt_err) {
addCOMPS <- function(mattop, tmat, doi_peaks, tmp, mz_err, rt_err) {

  ## for every component in the selected matching table
  for(cmp in na.omit(mattop$comp)) {

    ## take single DOI COMP
    component <- mattop %>%
      filter(comp == cmp)

    ####--- (A) check if this DOI COMP was assigned with a CID
    if (any(!is.na(component$topCOMP))) {
      update <- addGROUPED(component = component, tmat = tmat, tmp = tmp, mz_err = mz_err, rt_err = rt_err, doi_peaks = doi_peaks)

    } else {
      ####---- (B) if DOI component was not assigned with a CID, add its peaks as new
      update <- addUNGROUPED(component = component, tmat = tmat, tmp = tmp, mz_err = mz_err, rt_err = rt_err, doi_peaks = doi_peaks)

    }

    ####---- universal part for updating tmp
    ## extract from list
    utmp <- update$utmp
    doi_peaks <- update$doi_peaks

    ## save a copy of the template without the template peaks matched by mz/rt window
    ttmp <- tmp %>%
      filter(!pid %in% (utmp %>% pull(pid)))

    ## bind peaks to TMP
    tmp <- bind_rows(ttmp, utmp)

  }

  return(list("doi_peaks" = doi_peaks, "tmp" = tmp))

}


####---- If any of component's peaks were grouped to any template's CID, add all of its peaks to that CID
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

###---- if component was not assigned/grouped with any template's CID, add its peaks to template as new
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

  if( nrow(cmat) != nrow(utmp)) { stop(paste("un-assigned DOI component peaks are not added to TMP correctly! Check p:", p)) }

  ## assign CIDs for the checked DOI peaks
  doi_peaks[which(doi_peaks$pno %in% (cmat %>% filter(target == T) %>% pull(pno))),] <- utmp %>%
    filter(pno %in% (cmat %>% filter(target == T) %>% pull(pno))) %>%
    select(pno, pid, cid, clid)

  return(list("doi_peaks" = doi_peaks, "utmp" = utmp))

}

####---- Select peaks if TMP peak matches to multiple keys, decide which one is closer in MZ
selectPEAK <- function(t, cmat){

  if (any(t$target == F) & all(!is.na(t$key))) {

    closest <- t %>%
        filter(target == F) %>%
        group_by(pid) %>%
        mutate(key_to_match = key) %>%
        ## find mz difference between the KEY and the matched tmp peaks
        mutate(dif = abs((cmat %>% ungroup() %>% filter(target == T, key == key_to_match ) %>% select(mz) %>% pull()) - mz)) %>%
        group_by(key) %>%
        filter(dif == min(dif)) %>%
        select(-c(key_to_match, dif))

  } else {

    closest <- t
  }

  return(closest)


}

