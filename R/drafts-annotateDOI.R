####---- get peakgroup clusters
## requires columns 'peakgr', 'peakid' and rt/mz values
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


####---- find matches for a peak from one dataset in another
## use peak's and dataset's mz and rt regions only
matchPEAK <- function(p, dt) {
  mat <- dt %>%
    ## Using regions of both peaks. Regions defined by user mz and rt error ranges
    filter(between(mz_l, p$mz_l, p$mz_h) | between(mz_h, p$mz_l, p$mz_h)) %>%
    filter(between(rt_l, p$rt_l, p$rt_h) | between(rt_h, p$rt_l, p$rt_h)) %>%
    mutate(target_peakid = p$peakid, target_peakgr = p$peakgr) %>%
    mutate(target_peakgrcls = p$peakgrcls) %>%
    # select(target_peakid, target_peakgr, names(dt))
    select(target_peakid, target_peakgr, target_peakgrcls, peakid, peakgr, peakgrcls, mz, rt, into)

  return(mat)

}

####---- Prepare datafile-of-interest for grouping with the template
getDOI <- function(file = NULL) {
  if (is.null(file)) stop("file is required")
  if (!file.exists(file)) stop("incorrect filepath in the provided.\n")
  doi <- read.csv(file, stringsAsFactors = F)
  required_colnames <- c("peakid", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "peakgr")
  if (any(!required_colnames %in% colnames(doi))) stop("file is missing columns: ",
                                                       paste0(required_colnames[which(!required_colnames %in% colnames(dt))], collapse = ", "))

  doi <- getCLUSTS(dt = doi)
  return(doi)
}

####---- compare all matched template's peaks and find the top matches for every target peakgroup
getTOPmatches <- function(mat, target, tmp, bins){

  ## (1) add peaks, that belong to the same cluster as the matched cluster(s)
  ## bind rows, rather than extract from template, so that mat columns 'target_peakgr', 'target_peakid' would be retained
  allmat <- bind_rows(mat,
                      tmp %>%
                        filter(peakgrcls %in% (mat %>% distinct(peakgrcls) %>% pull(peakgrcls))) %>% ## add additional cluster peaks
                        filter(!peakid %in% (mat %>% distinct(peakid) %>% pull(peakid))) %>% ## but not those, that are already in the mat dataframe
                        select(peakid, peakgr, peakgrcls, mz, rt, into )) ## minimise dataframe to the most important columns

  ## check similarity of each target peakgroup against all its matches in the template
  matcos <- target %>%
    rename(target_peakgr = peakgr, target_peakgrcls = peakgrcls) %>% ## must have a different column names to be retained in the final dataframe output
    ## for each peakgroup in the target
    group_by(target_peakgr, target_peakgrcls) %>%
    do(compareSPECTRA(t = ., allmat = allmat, bins = bins))

  ## compare clusters,  extract top matches
  mattop <- compareCLUSTERS(matcos = matcos, allmat = allmat, target = target)

  return(mattop)

}

####---- compare spectra of the target peakgroup and all matching template peakgroups
compareSPECTRA <- function(t, allmat, bins) {

  ## extract matches for the single target peakgroup
  tmat <- allmat %>% filter(target_peakgr == unique(t$target_peakgr))

  if (nrow(tmat) > 0) {
    tmat <- allmat %>%
      filter(peakgrcls %in% unique(tmat$peakgrcls))

    ## build total peak table for both the target peaks and all template matches
    all <- bind_rows(t %>%
                       select(peakid, peakgr = target_peakgr, peakgrcls = target_peakgrcls, mz, into) %>%
                        mutate(dt = "target"),
                      tmat %>%
                        select(peakid, peakgr, peakgrcls, mz, into) %>%
                        mutate(dt = "matched"))

    ## generate mz bins using all peaks
    breaks <- data.frame(breaks = seq(from = min(all$mz), to = max(all$mz), by = bins)) %>%
      mutate(bin = row_number())
    all <- all %>%
      mutate(bin = findInterval(mz, breaks$breaks))

    ## build target vector of the target spectrum
    target_vec <- full_join(breaks,
                           all %>%
                              filter(dt == "target") %>%
                              group_by(bin) %>% filter(row_number() == 1) %>% ## remove duplicating duplicating bins from different PIDS
                              select(bin, into),
                            by = c("bin")) %>%
      mutate(into = ifelse(is.na(into), 0, into)) %>%
      mutate(into = scaleVEC(into)) %>% ## scale vector to unit length of 1
      pull(into)

    ## get the cosine of the angle between vectors of the target and the matching template peakgroups
    cos <- all %>%
      filter(dt == "matched") %>%
      group_by(peakgrcls, peakgr) %>%
      do(data.frame(cos = getCOS(breaks = breaks, target_vec = target_vec, matched = .))) %>% ## named argument cos would become a list if not put inside data.frame()
      ungroup()
  } else {
    cos <- data.frame(peakgrcls = NA, peakgr = NA, cos = NA)
  }

  return(cos)
}

####---- Cosine estimation
getCOS <- function(breaks, target_vec, matched){

  ## build vector
  matched_vec <- full_join(breaks,
                           matched %>%
                            group_by(bin) %>% filter(row_number() == 1) %>% ## remove duplicating duplicating bins from different PIDS
                            select(bin, into),
                          by = c("bin")) %>%
    mutate(into = ifelse(is.na(into), 0, into)) %>% ## if mz bin has an intensity, use it, otherwise set to 0
    mutate(into = scaleVEC(into)) %>% ## scale vector to unit length of 1
    pull(into)

  ## find the cosine of the angle
  cos_angle <- (sum(target_vec * matched_vec) )  / ( (sqrt(sum(target_vec * target_vec)))  * ( sqrt(sum(matched_vec * matched_vec)) ) )
  return(cos_angle)

}

####---- Scale vector to unit length
scaleVEC <- function(x) {
  x / ( sqrt(sum(x * x)) )
}

####---- Compare clusters
compareCLUSTERS <- function(matcos, allmat, target) {

 ## convert to ranks based on cosine
 mattop <- matcos %>%
   filter(cos > 0) %>%
   ungroup() %>%
   ## if any pairs have identical cosines (only possible when running with simulated datasets), mark for further evaluation
   add_count(cos) %>%
   mutate(duplicated = ifelse(any(n > 1), TRUE, FALSE)) %>%
   ## for each template and target peakgr separately find the order of matching pairs based on cosine
   group_by(peakgr) %>%
   mutate(template_rank = ifelse(duplicated == FALSE, row_number(desc(cos)), NA)) %>%
   group_by(target_peakgr) %>%
   mutate(target_rank = ifelse(duplicated == FALSE, row_number(desc(cos)), NA)) %>%
   ## mark pairs that are TOP
   group_by(peakgr) %>%
   mutate(top_rank =
            ## if all cosines are unique
            if (all(duplicated == FALSE)) {
              ## if any peakgr pair is the top by both template, and target aproach
              if (any(template_rank == 1 & target_rank == 1)) {
                ifelse(template_rank == 1 & target_rank == 1, TRUE, FALSE)
                } else {
                  NA
                }
              } else {
               ## further analysis
               "further"
             }) %>%
   group_by(target_peakgr) %>%
   mutate(top =
            ## if all cosines are unique
            if (all(duplicated == FALSE)) {
              ## (A) if TOP was NOT assigned to any peakgr for this target_peakgr, take next best that is not assigned to anything yet
              if (any(top_rank[!is.na(top_rank)] == FALSE)) {
                ifelse(is.na(top_rank) & template_rank == min(template_rank[is.na(top_rank)]), TRUE, FALSE)
                } else {
                  ## (B) if TOP was assigned to any peakgr for this target_peakgr, mark it
                  if (any(top_rank[!is.na(top_rank)] == TRUE)) {
                    ifelse(top_rank == TRUE, TRUE, FALSE)
                    } else {
                      NA
                    }
                  }
              } else {
                ## further investigation
                "further"
                })

 if (any(na.omit(mattop$top == "further"))) { warning("identical cosines were found!") }

 ## for each template peakgroup, retain only one TOP target peakgroup (if any)
 mattop <- mattop %>%
   group_by(peakgr) %>%
   filter(top == T) %>%
   ungroup() %>%
   select(target_peakgr, target_peakgrcls, peakgr, peakgrcls, cos)

 ## add target peakgroups that did not have a match, will be needed when updating doi dataframe
 mattop <- matcos %>%
   group_by(target_peakgr, peakgr) %>%
   mutate(top = ifelse(target_peakgr %in% mattop$target_peakgr & peakgr %in% mattop$peakgr, TRUE, FALSE))

  return(mattop)

}

####---- add DOI peaks to TMP
## function returns updated tmp and doi dataframes
addPEAKS <- function(mattop, tmp) {

  ## for every target peakgr
  for(pkg in unique(na.omit(mattop$target_peakgr))) {

    peakgroup <- mattop %>%
      filter(target_peakgr == pkg)

    ####--- (A) if this target peakgroup was assigned (top == TRUE) and also had a higher cos than any previously assigned DOI peakgroup
    if (any(!is.na(peakgroup$top))) {

      ## check if assigned peakgroup was not grouped with others from the same DOI
      previous <- tmp %>%
        filter(peakgr == peakgroup$peakgr) %>%
        select(peakid, peakgr, contains("cos"), contains("target_peakid")) %>% ## use contains() if tmp doesn't have column 'cos' yet (first grouping)
        mutate(cos = ifelse(ncol(.) > 2 , cos, 0)) %>% ## if cid doesn't have cos (first grouping), assign cos with 0
        mutate(cos = ifelse(any(!is.na(cos)), na.omit(cos), 0))

      peakgrouptop <- ifelse( peakgroup %>% filter(top == T) %>%  pull(cos) > (unique(previous$cos)), T, F)

      if (peakgrouptop) {
        update <- addGROUPED(component = component, target = target, allmat = allmat, tmpo = tmpo, tmp = tmp, doi = doi, previouscomp = previouscomp, mz_err = mz_err, rt_err = rt_err, doi_peaks = doi_peaks)

      } else {
        ####---- (B) otherwise, add its peaks as new
        update <- addUNGROUPED(component = component, target = target, allmat = allmat, tmp = tmp, mz_err = mz_err, rt_err = rt_err, doi_peaks = doi_peaks)

      }

    } else {
      ####---- (B) otherwise, add its peaks as new
      update <- addUNGROUPED(component = component, target = target, allmat = allmat,  tmp = tmp, mz_err = mz_err, rt_err = rt_err, doi_peaks = doi_peaks)

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

    ## remove pids marked in 'pid_to_remove'
    tmp <- tmp %>%
      filter(pid_to_remove == F | is.na(pid_to_remove))

  }

  return(list("doi_peaks" = doi_peaks, "tmp" = tmp))

}

####---- If any of component's peaks were grouped to any template's CID, add all of its peaks to that CID
addGROUPED <- function(peakgroup, target, allmat, tmpo, tmp, doi, previous, mz_err, rt_err, doi_peaks) {

  ## take the TMP peakgroup, which was assigned to the component
  peakgroup <- peakgroup %>%
    filter(top == T)

  pkg <- target %>%
    filter(peakgr == peakgroup$target_peakgr)

  ####---- (A) if this CID was previously grouped with other COMP, and now has to be over-written by current COMP
  # if (unique(previous$cos) > 0) {
  #
  #   ####---- (A1) remove peaks in TMP coming from previous COMP
  #   ### re-assign previous COMP to new CID and PIDs
  #   ### here mz/rt/into come from DOI and can be readily merged to tmp
  #   ### assign pids T or F in column 'pid_to_remove' to mark for row removal from tmp once merged with tmp
  #   ### get a new CID for previous COMP
  #   cid_id <- max(tmp$cid) + 1
  #
  #   ## table with old pids, which are marked 'pid_to_remove' = T for removal once merged with tmp
  #   ptmp <- previouscomp %>%
  #     filter(!pid %in% tmpo$pid) %>%
  #     mutate(pid_to_remove = TRUE)
  #
  #   ## duplicate table with new pids (and all grouping metadata), which will be kept in the tmp
  #   ptmp <- previouscomp %>%
  #     filter(pno %in% doi$pno) %>%
  #     group_by(pno) %>%
  #     mutate(pid_to_remove = F) %>%
  #     ## for peaks coming from DOI, assign original DOI mz/rt values
  #     mutate(mz = ifelse(pno %in% doi$pno, doi[which(doi$pno == pno),"mz"], NA)) %>%
  #     mutate(rt = ifelse(pno %in% doi$pno, doi[which(doi$pno == pno),"rt"], NA)) %>%
  #     mutate(into = ifelse(pno %in% doi$pno, doi[which(doi$pno == pno),"into"], NA)) %>%
  #     mutate(pid = NA, cid = cid_id, cos = NA) %>%
  #     mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
  #     ungroup() %>%
  #     bind_rows(., ptmp)
  #
  #   ## change CIDs for peaks from the previous COMP
  #   doi_peaks <- right_join(doi_peaks %>%
  #                             filter(pno %in% (previouscomp %>% filter(!is.na(pno)) %>% pull(pno))) %>%
  #                             mutate(cid = cid_id),
  #                           doi_peaks, by = "pno") %>%
  #     mutate(cid = ifelse(is.na(cid.x), cid.y, cid.x)) %>%
  #     select(pno, cid)
  #
  #
  #   ####---- (A2) update template with current COMP
  #   ## take TMP peaks that were returned by matched CLID
  #   matcomp <- allmat %>%
  #     filter(
  #       (comp == component$comp & cid == component$cid) |
  #         (is.na(comp) & cid == component$cid)) ## peaks not matched my mz/rt, but by CID
  #
  #   ## merge TMP peaks that match to the same PNO
  #   ## for each PNO that matches multiple tmp peaks, find the closest peak
  #   matcomp <- matcomp %>%
  #     group_by(pno) %>%
  #     do(selectPEAK(m = ., target = target)) %>%
  #     ungroup()
  #
  #   ## select and average can be put into one function!! future develop
  #   ## for each PNO, average the mz/rt across all matches. Returns all target peaks and their matches (if any)
  #   ntmp <- targetcomp %>%
  #     group_by(pno) %>%
  #     do(getMEDIAN(t = ., m = matcomp)) %>%
  #     ungroup()
  #
  #   ntmp <- ntmp %>%
  #     bind_rows(.,
  #               ## update TMP peaks that were not matched
  #               matcomp %>%
  #                 filter(!pid %in% ntmp$pid) %>%
  #                 mutate(tmp = TRUE, comp = component$comp, cls = unique(targetcomp$cls))) %>%
  #     mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
  #     mutate(cid = component$cid, clid = component$clid,
  #            cos = unlist(component$cos), topCOMP = component$topCOMP, tmpCLS_order = component$tmpCLS_order, targetCLS_sc = component$targetCLS_sc, tmpCLS_sc = component$tmpCLS_sc) %>%
  #     mutate(pid_to_remove = FALSE)
  #
  #   utmp <- bind_rows(ntmp, ptmp)
  #
  # } else {

    ####---- (B) if this peakgr was NOT previously grouped, simply add/merge peaks to it
    matpkg <- allmat %>%
      filter(
        (peakgr == peakgroup$peakgr &target_peakgr == peakgroup$target_peakgr) |
          (peakgr == peakgroup$peakgr & is.na(target_peakgr))) ## peaks not matched my mz/rt, but by peakgr

    ## STOPPED HERE
    ## merge TMP peaks that match to the same PNO
    ## for each PNO that matches multiple tmp peaks, find the closest peak
    matpkg <- matpkg %>%
      group_by(peakid) %>%
      do(selectPEAK(m = ., target = target)) %>%
      ungroup()

    ## select and average can be put into one function!! future develop

    ## for each PNO, average the mz/rt across all matches. Returns all target peaks and their matches (if any)
    utmp <- targetcomp %>%
      group_by(pno) %>%
      do(getMEDIAN(t = ., m = matcomp)) %>%
      ungroup()

    ## add matching metadata
    utmp <- utmp %>%
      bind_rows(.,
                ## update TMP peaks that were not matched
                matcomp %>%
                  filter(!pid %in% utmp$pid) %>%
                  mutate(tmp = TRUE, comp = component$comp, cls = unique(targetcomp$cls))) %>%
      mutate( mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
      mutate(cid = component$cid, clid = component$clid,
             cos = unlist(component$cos), topCOMP = component$topCOMP, tmpCLS_order = component$tmpCLS_order, targetCLS_sc = component$targetCLS_sc, tmpCLS_sc = component$tmpCLS_sc)
  }


  ## give PIDs for newly added/re-added peaks
  if (nrow(utmp %>% filter(is.na(pid))) > 0) {
    pid_ids <- (max(tmp$pid) + 1):( max(tmp$pid) + nrow(utmp %>% filter(is.na(pid))) )
    utmp[which(is.na(utmp$pid)), "pid"] <- pid_ids
  }

  ####----(4) assign CIDs for the current COMP
  doi_peaks <- right_join(utmp %>%
                            filter(pno %in% (targetcomp %>% pull(pno))) %>%
                            filter(cid == component$cid) %>%
                            select(pno, cid),
                          doi_peaks, by = "pno") %>%
    mutate(cid = ifelse(is.na(cid.x), cid.y, cid.x)) %>%
    select(pno, cid)

  return(list("doi_peaks" = doi_peaks, "utmp" = utmp))

}






######## major function ########
####---- Add DOI peaks to the DB annotations table
## add DATAFILE (doi) to the the TEMPLATE (tmp)
## template could be either the master peak table, or a DATABASE (if building template for the first time)
## returns a template
addDOI <- function(tmp, doi, mz_err = 0.01, rt_err = 2, bins = 0.1) {

  ## (A) prepare doi for matching against template/db
  doi_full <- getDOI(file = doi)
  doi <- doi_full %>%
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
    mutate(chemid = NA, dbid = NA, added = F) %>%
    select(peakid, mz, mz_l, mz_h, rt, rt_l, rt_h, into, peakgr, peakgrcls, chemid, dbid, added)

  ## (B) prepare template for grouping with the doi
  # db is the template where doi peaks will be added to
  tmp <- getCLUSTS(dt = tmp)
  tmp <- tmp %>%
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err)

  ## for every doi peak, find matching tmp peaks (if any)
  while (any(doi$added == F)) {

    ## get the first peak among the un-annotated peaks
    p <- doi %>% filter(added == F) %>% slice(1) %>% pull(peakid)
    target <- doi %>% filter(peakid == p)
    target <- doi %>% filter(peakgrcls == target$peakgrcls)

    ## matching by mz/er window
    mat <- target %>%
      group_by(peakid) %>%
      do(matchPEAK(p = ., dt = tmp)) %>%
      ungroup()

    ## extract top matches for every target peakgr
    mattop <- getTOPmatches(mat = mat, target = target, tmp = tmp)

    ## add annotated/unmatched DOI peaks to template
    update <- addPEAKS()

    }



}
