####---- get peakgroup clusters based on rt values
## requires columns 'peakgr','peakid','rt'
## returns original dt with one new column 'peakgrcls'
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
      cluster <- dplyr::bind_rows(cluster,
                                  dt %>%
                                    filter(peakgr %in% (cluster %>% distinct(peakgr) %>% pull(peakgr))) %>%
                                    filter(!peakid  %in% (cluster %>% distinct(peakid) %>% pull(peakid))))

      ## update cluster ID assignment table
      cls[which(cls$peakgr %in% unique(cluster$peakgr)),"peakgrcls"] <- rep(id, length(unique(cluster$peakgr)))

    } else { next }

  }

  dt <- dplyr::full_join(dt, cls, by = c("peakgr"))

  ## order peak-groups by their peakgrcomplexity
  ## and peaks intensity
  dt_c <- dt %>%
    group_by(peakgr) %>%
    summarise(n = n()) %>%
    dplyr::arrange(desc(n)) %>%
    ungroup()

  dt <- dt %>%
    dplyr::arrange(factor(peakgr, levels = dt_c$peakgr), peakgr, peakid)

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
    select(target_peakid, target_peakgr, names(dt))
    # select(target_peakid, target_peakgr, target_peakgrcls, peakid, peakgr, peakgrcls, mz, rt, into)

  return(mat)

}

####---- Prepare datafile-of-interest for alignmemt
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
   select(target_peakgr, target_peakgrcls, peakgr, peakgrcls, cos, top)

 ## add target peakgroups that did not have a match, will be needed when updating doi dataframe
 mattop <- mattop %>%
   bind_rows(.,
             matcos %>%
               filter(!target_peakgr %in% mattop$target_peakgr) %>%
               mutate(peakgr = NA, peakgrcls = NA, cos = NA, top = FALSE))

  return(mattop)

}

####---- add DOI peaks to TMP
## function returns updated tmp and doi dataframes
addPEAKS <- function(mattop, mat, target, tmp, itmp, doi) {

  ## for every target peakgr
  for(pkg in unique(na.omit(mattop$target_peakgr))) {

    matpkg <- mattop %>%
      filter(target_peakgr == pkg, top == TRUE)

    ####--- (A) if this target peakgr was assigned with a template peakgr
    if (nrow(matpkg) > 0) {
      update <- addGROUPED(matpkg = matpkg, target = target, mat = mat, itmp = itmp, doi = doi)

    } else {
      ####---- (B) otherwise, add its peaks as new
      update <- addUNGROUPED()
    }

    ####---- universal part for updating tmp
    utmp <- update$utmp

    ## save a copy of the tmp without the tmp peaks matched by mz/rt window
    ttmp <- itmp %>%
      filter(!peakid %in% (utmp %>% pull(peakid)))

    ## bind peaks to tmp
    itmp <- bind_rows(ttmp, utmp)

  }

  return(list("doi" =  update$doi, "itmp" = itmp))

}




####---- If any of this target peakgr are grouped to template peakgr, add all of its peaks to that peakgr
addGROUPED <- function(matpkg, target, mat, itmp, doi) {

  ## extract target doi peaks that must be added to tmp
  doipeaks <- target %>%
    filter(peakgr == matpkg$target_peakgr)

  ## extract tmp peaks that must be merged/averaged with doi peaks
  tmppeaks <- mat %>%
    filter(peakgr ==  matpkg$peakgr,
           target_peakgr == matpkg$target_peakgr)

  ## extract template peakgr
  previous <- itmp %>%
    filter(peakgr == matpkg$peakgr) %>%
    mutate(cos = ifelse(!is.na(cos), cos, 0))

  ####---- (A) if previously grouped doi peaks were a better fit than current doi peaks
  if (unique(previous$cos) > matpkg$cos) {

    ## retain previously grouped doi peaks in the tmp unmodified
    ## add current doi peaks to a newly generated peakgr instead
    peakid_new <- (max(itmp$peakid) + 1):(max(itmp$peakid) + nrow(doipeaks))
    peakgr_new <- max(itmp$peakgr) + 1

    utmp <- doipeaks %>%
      mutate(doi_peakid = peakid, doi_peakgr = peakgr, doi_peakgrcls = peakgrcls, cos = NA) %>%
      mutate(peakid = peakid_new, peakgr = peakgr_new, peakgrcls = NA) %>%
      select(peakid, mz, rt, into, peakgr, peakgrcls, doi_peakid, doi_peakgr, doi_peakgrcls, cos)

    } else {

      ####---- (A2) if previously grouped doi peaks were a worse fit and have to be removed from peakgr
      if (unique(previous$cos) > 0) {

        ## update previously added doi peaks
        previous_peakid <- previous %>%
          filter(!is.na(doi_peakid)) %>%
          pull(doi_peakid)
        peakid_new <- (max(itmp$peakid) + 1):(max(itmp$peakid) + length(previous_peakid))
        peakgr_new <- max(itmp$peakgr) + 1
        previouspeaks <- doi %>% ## use original doi values, instead of merged itmp values
          filter(peakid %in% previous_peakid) %>%
          mutate(doi_peakid = peakid, doi_peakgr = peakgr, doi_peakgrcls = peakgrcls) %>%
          mutate(peakid = peakid_new, peakgr = peakgr_new, peakgrcls = NA, cos = NA) %>%
          select(peakid, mz, rt, into, peakgr, peakgrcls, doi_peakid, doi_peakgr, doi_peakgrcls, cos)
      } else {
        ## prepare tmp peaks that are part of the same tmp peakgr
        previouspeaks <- data.frame(stringsAsFactors = F)
      }

      ##  merge current doi peaks with the selected tmp peakgr
      ## for every doi peak:
      ## find the closest peak in the template (if multiple matches)
      ## average mz and rt across tmp and doi
      newpeaks <- doipeaks %>%
        group_by(peakid) %>%
        do(averagePEAKS(p = ., tmppeaks = tmppeaks, matpkg = matpkg)) %>%
        ungroup() %>%
        ## add chemical db metadata
        mutate(chemid = unique(tmppeaks$chemid), dbid = unique(tmppeaks$dbid), dbname = unique(tmppeaks$dbname))
      peakid_new <- (max(itmp$peakid) + 1):(max(itmp$peakid) + nrow(newpeaks))
      newpeaks <- newpeaks %>%
        mutate(peakid = ifelse(is.na(peakid), peakid_new, peakid))

      utmp <- bind_rows(previouspeaks, newpeaks) %>%
        bind_rows(previous %>%
                    filter(!peakid %in% newpeaks$peakid) %>%
                    mutate(doi_peakgr = unique(newpeaks$doi_peakgr),
                           doi_peakgrcls = unique(newpeaks$doi_peakgrcls),
                           cos = unique(newpeaks$cos)))

  }

  doi[which(doi$peakid %in% na.omit(utmp$doi_peakid)),
      c("tmp_peakid", "tmp_peakgr", "chemid", "dbid", "dbname", "cos")] <- utmp %>%
    filter(!is.na(doi_peakid)) %>%
    select(peakid, peakgr, chemid, dbid, dbname, cos)

  doi[which(doi$peakid %in% na.omit(utmp$doi_peakid)),"added"] <- rep(TRUE, length(na.omit(utmp$doi_peakid)))

  return(list("doi" = doi, "utmp" = utmp))

}












####---- Select peaks if TMP peak matches to multiple target peaks, decide which one is closer in MZ
## iterated over every doi peakid (p)
averagePEAKS <- function(p, tmppeaks, matpkg){

  ## select closest tmp peak
  closest <- tmppeaks %>%
    filter(target_peakid == p$peakid) %>%
    group_by(peakid) %>%
    ## find mz difference between the target peakid and the matched tmp peakid
    mutate(dif = abs(p$mz - mz)) %>%
    filter(dif == min(dif)) %>%
    ## if more than 1 peakid has the same MZ difference to target, take the most intense peakid
    filter(into == max(into))

  if (nrow(closest) > 0) {

    if (nrow(closest) > 1) stop("averagePEAKS fails to select a single peak")

    ## change tmp peak's mz,rt to average, into to the doi
    peaks <- closest %>%
      mutate(mz = median(c(mz, p$mz)), rt = median(c(rt, p$rt)), into = p$into) %>%
      mutate(doi_peakid = p$peakid, doi_peakgr = p$peakgr, doi_peakgrcls = p$peakgrcls, cos = matpkg$cos)

  } else {
    ## if target PNO doesn't have a match, add it as a new peak to the template
    peaks <- p %>%
      mutate(doi_peakid = peakid, doi_peakgr = peakgr, doi_peakgrcls = peakgrcls, cos = matpkg$cos) %>%
      mutate(peakid = NA, peakgr = matpkg$peakgr, peakgrcls = matpkg$peakgrcls) ## peakid will be generated later on during tmp merging
  }
  peaks <- peaks %>%
    select(peakid, mz, rt, into, peakgr, peakgrcls, doi_peakid, doi_peakgr, doi_peakgrcls, cos)

  return(peaks)
}



