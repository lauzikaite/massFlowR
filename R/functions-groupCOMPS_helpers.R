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

  if(nrow(mat) == 0 ) {

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

  # multi <- tmat %>%
  #   filter(target == F) %>%
  #   group_by(comp) %>%
  #   summarise(target_comp = paste0(unique(cid), collapse = ","),
  #             target_comp_n = length(unique(cid)))


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

    # multi <- tmat %>%
    #   filter(target == F) %>%
    #   group_by(pid) %>%
    #   summarise(target_comp = paste0(unique(comp), collapse = ","),
    #             target_key = paste0(unique(key), collapse = ","),
    #             target_key_n = length(unique(key)))

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
                        group_by(cid, cls) %>%
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
    group_by(cls, cid) %>%
    summarise(rtmed = median(rt)) %>%
    arrange(rtmed) %>%
    group_by(cls) %>%
    mutate(order = row_number()) %>%
    ungroup()

  ## add components order in the cluster
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
    group_by(cls) %>%
    arrange(target_order) %>%
    filter(topCOMP == T) %>%
    summarise(tmpCLS_order = all(tmp_order == cummax(tmp_order)),
              targetCLS_sc = sum(!is.na(topCOMP))/ max(o_target$order),
              tmpCLS_sc =  sum(!is.na(topCOMP))/ max(tmp_order))
  ## add summaries
  mattop <- full_join(mattop %>% select(comp, cls, cid, cos, topCOMP),
                      matcls, by = c("cls"))

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
  update <- updateTMP(mattop = mattop, tmat = tmat, doi_peaks = doi_peaks, tmp = tmp, mz_err = mz_err, rt_err = rt_err)

  return(list("doi_peaks" = update$doi_peaks, "tmp" = update$tmp))

}









####---- Update template with selected matches
updateTMP <- function(mattop, tmat, doi_peaks, tmp, mz_err, rt_err) {

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







####---- Update template by adding new peaks that were not present in template before
# addUNGROUPEDnew <- function(cmat, component, tmp, mz_err, rt_err) {
#
#   ## assign new CID and PID for the unmatched target component
#   cid_id <- max(tmp$cid) + 1
#   pid_id <- (max(tmp$pid) + 1):(max(tmp$pid) + nrow(cmat))
#
#   utmp <- cmat %>%
#     mutate(pid = pid_id,
#            clid = NA, # temporalily assign cluster ID with NA, since ids will be generated again at the end of tmp grouping round
#            mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err,
#            tmp = FALSE,
#            ## since keys are not pasted, numeric format would conflict with other scenarios
#            key = as.character(key), key_max = as.character(key_max),
#            ## DOI's info
#            comp = cid,
#            cos = NA,
#            ## TMP's info comes last since this is erasing
#            cid = cid_id,
#            tmpCLS_order = NA,
#            targetCLS_sc = component$targetCLS_sc,
#            tmpCLS_sc = NA) %>%
#     select(pid, mz, rt, into, cid, clid, mz_l, mz_h, rt_l, rt_h, tmp, pno, key, key_max, comp, cls, cos, tmpCLS_order, targetCLS_sc, tmpCLS_sc)
#
#   return(utmp)
# }

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

