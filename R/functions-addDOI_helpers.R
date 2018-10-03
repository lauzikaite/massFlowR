####---- Check and load file for alignment
checkFILE <- function(file = NULL) {
  if (!file.exists(file)) stop("incorrect filepath for: ", file)
  doi <- read.csv(file, stringsAsFactors = F)
  required_colnames <- c("peakid", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "peakgr")
  if (any(!required_colnames %in% colnames(doi))) {
    stop("file ", file, "is missing columns: \n",
         paste0(required_colnames[which(!required_colnames %in% colnames(dt))], collapse = ", "))
  }
  return(doi)
}


####---- add columns to dataframe, fill with NA
addCOLS <- function(dt, cname) {
  add <- cname[!cname %in% names(dt)]
  if (length(add) != 0) dt[add] <- NA
  return(dt)
}

####---- get peakgroup clusters based on rt values
## requires columns 'peakgr','peakid','rt'
## returns original dt with one new column 'peakgrcls'
getCLUSTS <- function(dt){

  ## data frame to store assigned peak-group-clusters' ids
  cls <- data.frame(peakgr = unique(dt$peakgr), peakgrcls = NA, stringsAsFactors = F)

  ## for every peak group
  for (pg in unique(dt$peakgr)) {

    ## if peak-group was not assigned to a cluster yet
    if (is.na(cls[which(cls$peakgr == pg),"peakgrcls"])) {

      ## get unique peak-group-cluster id
      id <- ifelse(all(is.na(cls$peakgrcls)), 1, max(na.omit(cls$peakgrcls)) + 1)

      ## find all other peaks in the same RT region as the peak-group
      peaks <- dt %>%
        filter(.data$peakgr == pg)
      cluster <- dt %>%
        filter(.data$peakgr %in% (cls %>% filter(is.na(.data$peakgrcls)) %>% select(.data$peakgr) %>% pull())) %>% # only search between peak groups that were not clustered yet
        filter(between(.data$rt, min(peaks$rt), max(peaks$rt))) # use RT region defined by the min and max values

      ## add other peaks in the extracted peak-group
      cluster <- dplyr::bind_rows(cluster,
                                  dt %>%
                                    filter(.data$peakgr %in% (cluster %>% distinct(.data$peakgr) %>% pull(.data$peakgr))) %>%
                                    filter(!.data$peakid  %in% (cluster %>% distinct(.data$peakid) %>% pull(.data$peakid))))

      ## update cluster ID assignment table
      cls[which(cls$peakgr %in% unique(cluster$peakgr)),"peakgrcls"] <- rep(id, length(unique(cluster$peakgr)))

    }
    else next

  }

  dt <- dplyr::full_join(dt, cls, by = c("peakgr"))

  ## order peak-groups by their peakgrcomplexity
  ## and peaks intensity
  dt_c <- dt %>%
    group_by(.data$peakgr) %>%
    summarise(n = n()) %>%
    arrange(desc(.data$n)) %>%
    ungroup()

  dt <- dt %>%
    dplyr::arrange(factor(.data$peakgr, levels = dt_c$peakgr), .data$peakgr, .data$peakid)

  return(dt)

}


####---- find matches for a peak from one dataset in another
## use peak's and dataset's mz and rt regions only
matchPEAK <- function(p, dt) {
  mat <- dt %>%
    ## Using regions of both peaks. Regions defined by user mz and rt error ranges
    filter(between(.data$mz_l, p$mz_l, p$mz_h) | between(.data$mz_h, p$mz_l, p$mz_h)) %>%
    filter(between(.data$rt_l, p$rt_l, p$rt_h) | between(.data$rt_h, p$rt_l, p$rt_h)) %>%
    mutate(target_peakid = p$peakid, target_peakgr = p$peakgr) %>%
    mutate(target_peakgrcls = p$peakgrcls) %>%
    select(.data$target_peakid, .data$target_peakgr, names(dt))
  return(mat)
}

####---- compare all matched template's peaks and find the top matches for every target peakgroup
getTOPmatches <- function(mat, target, tmp, bins, add_db, db_thrs) {

  if (nrow(mat) > 0) {
    ## (1) add peaks, that belong to the same cluster as the matched cluster(s)
    ## bind rows, rather than extract from template, so that mat columns 'target_peakgr', 'target_peakid' would be retained
    allmat <- bind_rows(mat,
                        tmp %>%
                          filter(.data$peakgrcls %in% (mat %>% distinct(.data$peakgrcls) %>% pull(.data$peakgrcls))) %>% ## add additional cluster peaks
                          filter(!.data$peakid %in% (mat %>% distinct(.data$peakid) %>% pull(.data$peakid))) %>% ## but not those, that are already in the mat dataframe
                          select(.data$peakid, .data$peakgr, .data$peakgrcls, .data$mz, .data$rt, .data$into, .data$chemid)) ## minimise dataframe to the most important columns

    ## check similarity of each target peakgroup against all its matches in the template
    matcos <- target %>%
      rename(target_peakgr = .data$peakgr, target_peakgrcls = .data$peakgrcls) %>% ## must have a different column names to be retained in the final dataframe output
      ## for each peakgroup in the target
      group_by(.data$target_peakgr, .data$target_peakgrcls) %>%
      do(compareSPECTRA(t = ., allmat = allmat, bins = bins))

    ## compare all matching peakgrs in their clusters, extract top matches
    mattop <- compareCLUSTERS(matcos = matcos, add_db = add_db, db_thrs = db_thrs)

    } else {
      mattop <- data.frame(target_peakgr = unique(target$peakgr),
                           target_peakgrcls = unique(target$peakgrcls),
                           stringsAsFactors = F) %>%
        mutate(peakgrcls = NA, peakgr = NA, cos = NA, top = FALSE)
    }

  return(mattop)

}

####---- compare spectra of the target peakgroup and all matching template peakgroups
compareSPECTRA <- function(t, allmat, bins) {

  ## extract matches for the single target peakgroup
  tmat <- allmat %>% filter(target_peakgr == unique(t$target_peakgr))

  if (nrow(tmat) > 0) {
    tmat <- allmat %>%
      filter(.data$peakgrcls %in% unique(tmat$peakgrcls))

    ## build total peak table for both the target peaks and all template matches
    all <- bind_rows(t %>%
                       select(.data$peakid, peakgr = .data$target_peakgr, peakgrcls = .data$target_peakgrcls, .data$mz, .data$into) %>%
                        mutate(chemid = NA, dt = "target"),
                      tmat %>%
                        select(.data$peakid, .data$peakgr, .data$peakgrcls, .data$mz, .data$into, .data$chemid) %>%
                        mutate(dt = "tmp"))

    ## generate mz bins using all peaks
    breaks <- data.frame(breaks = seq(from = min(all$mz), to = max(all$mz), by = bins)) %>%
      mutate(bin = row_number())
    all <- all %>%
      mutate(bin = findInterval(.data$mz, breaks$breaks))

    ## build target vector of the target spectrum
    target_vec <- full_join(breaks,
                           all %>%
                              filter(dt == "target") %>%
                              group_by(.data$bin) %>% filter(row_number() == 1) %>% ## remove duplicating duplicating bins from different PIDS
                              select(.data$bin, .data$into),
                            by = c("bin")) %>%
      mutate(into = ifelse(is.na(.data$into), 0, .data$into)) %>%
      mutate(into = scaleVEC(.data$into)) %>% ## scale vector to unit length of 1
      pull(.data$into)

    ## get the cosine of the angle between vectors of the target and the matching template peakgroups
    cos <- all %>%
      filter(dt == "tmp") %>%
      group_by(.data$peakgrcls, .data$peakgr, .data$chemid) %>%
      do(data.frame(cos = getCOS(breaks = breaks, target_vec = target_vec, matched = .data))) %>% ## named argument cos would become a list if not put inside data.frame()
      ungroup()
  }
  else {
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
compareCLUSTERS <- function(matcos, add_db, db_thrs) {

  ## if aligning 1st doi in the study with the db, use similarity threshold (db_thrs)
  if (add_db == T) {
    mattop <- matcos %>%
      group_by(peakgr) %>%
      mutate(db = ifelse(!is.na(chemid), T, F)) %>%
      filter((db == T & cos > db_thrs) |
               (db == F & cos > 0)) %>%
      ungroup()
  } else {
    mattop <- matcos %>%
      filter(cos > 0) %>%
      ungroup()
  }

  if (nrow(mattop) > 0 ) {

   ## convert to ranks based on cosine
   mattop <- mattop %>%
     ## if pairs with the highest cosine have identical cosines (only possible when running with simulated datasets, to be investigated)
     add_count(cos) %>%
     mutate(duplicated = ifelse(
       n > 1 & cos == max(cos), T, F))
   if (any(mattop$duplicated == "TRUE")) {
     stop("identical cosines were found!")
   }

   ## (A) extract top pairs - each peakgr in the pair are ranked as 1 to each other
   mattop <- mattop %>%
     ## for each template and target peakgr separately, order the matches based on cosine
     group_by(peakgr) %>%
     mutate(template_rank = row_number(desc(cos))) %>%
     group_by(target_peakgr) %>%
     mutate(target_rank = row_number(desc(cos))) %>%
     ## for each target peakgr, take the top match (if available)
     group_by(target_peakgr) %>%
     mutate(top_rank =
              if (any(template_rank == 1 & target_rank == 1)) {
                ifelse(template_rank == 1 & target_rank == 1, TRUE, FALSE)
              } else {
                ## NA indicates that no TOP was assigned for this peakgr
                NA
              }) %>%
     ## for each tmp peakgr, mark if top was assigned, if not, mark top = NA
     group_by(peakgr) %>%
     mutate(top = if (any(na.omit(top_rank) == TRUE)) {
       ifelse(!is.na(top_rank), ifelse(top_rank == TRUE, TRUE, FALSE), FALSE)
       } else {
         NA
       })

   ## (B) extract next best pairs, by taking the pair with the highest cosine among the remaining tmp peakgr
   while (any(is.na(mattop$top))) {
     mattop <- maximiseASSIGNMENT(dt = mattop)
   }

   mattop <- mattop %>%
     group_by(target_peakgr) %>%
     filter(top == TRUE) %>%
     select(names(matcos), top) %>%
     ungroup()

 }

 ## add target peakgroups that did not have a match, will be needed when updating doi dataframe
 mattop <- mattop %>%
   bind_rows(.,
             matcos %>%
               filter(!target_peakgr %in% mattop$target_peakgr) %>%
               mutate(top = FALSE))

  return(mattop)

}


maximiseASSIGNMENT <- function(dt) {

  assigned <- dt %>%
    group_by(target_peakgr) %>%
    filter(!is.na(top))

  unassigned <- dt %>%
    group_by(peakgr) %>%
    filter(is.na(top)) %>%
    filter(!peakgr %in% assigned$peakgr[assigned$top == TRUE]) %>%
    filter(!target_peakgr %in% assigned$target_peakgr[assigned$top == TRUE])

  if (nrow(unassigned) > 0) {

  assigned_target_peakgrs <- assigned$target_peakgr
  assigned_target_peakgrs <- ifelse(length(assigned_target_peakgrs) == 0, 0, assigned_target_peakgrs)
  unassigned_cos <- unassigned$cos

  newly_assigned <- unassigned %>%
    group_by(peakgr) %>%
    mutate(top =
        ## if this peakgr has the highest cosine among the remaining peakgrs
        if (any(cos == max(unassigned_cos))) {
          ## if its target_peakgr hasn't been assigned already
          ifelse(cos == max(unassigned_cos), TRUE, FALSE)
        } else {
          NA
        })

  out <- bind_rows(assigned, newly_assigned)

  } else {
    out <- assigned
  }
  return(out)
}




####---- add DOI peaks to TMP
## function returns updated tmp and doi dataframes
addPEAKS <- function(mattop, mat, target, tmp, itmp, doi) {

  ## for every target peakgr
  for(pkg in unique(mattop$target_peakgr)) {

    matpkg <- mattop %>%
      filter(target_peakgr == pkg)

    ## take the single tmp top peakgr
    top_peakgr <- matpkg %>%
      filter(top == T) %>%
      pull(peakgr)

    ## take tmp peaks with the assigned peakgr
    ## these include both original tmp peaks, and doi peaks previously added to the tmp
    previous <- itmp %>%
      filter(peakgr == ifelse(length(top_peakgr) > 0, top_peakgr, 0))

    if (nrow(previous) > 0) {
      previous <- previous %>%
        mutate(cos = ifelse(!is.na(cos), cos, 0))
    } else {
      previous <- bind_rows(previous, data.frame(cos = 0, stringsAsFactors = F))
    }

    ####---- (A) if current assignment is better than previous
    if (any(matpkg$top == TRUE) & any(matpkg$cos > unique(previous$cos))) {
      update <- addGROUPED(matpkg = matpkg, target = target, mat = mat, previous = previous, itmp = itmp, doi = doi)
    } else {
      ####---- (B) if current assignment is worse than previous, or peakgr was not assigned, add these peaks as new peakgr in the tmp
      update <- addUNGROUPED(target = target, pkg = pkg, itmp = itmp, doi = doi)
    }

    ####---- universal part for updating tmp and doi
    doi <- update$doi
    utmp <- update$utmp
    ## bind new/modified peaks to itmp
    ttmp <- itmp %>%
      filter(!peakid %in% (utmp %>% pull(peakid))) %>%  ## save a copy of the tmp without the tmp peaks matched by mz/rt window
      filter(!doi_peakid %in% (utmp %>% filter(!is.na(doi_peakid)) %>% pull(doi_peakid))) ## save a copy of tmp without old doi peaks that now are updated
    itmp <- bind_rows(ttmp, utmp)
    
  }

  return(list("doi" = doi, "itmp" = itmp))

}




addGROUPED <- function(pkg, matpkg, target, mat, previous, itmp, doi) {

  matpkg <- matpkg %>%
    filter(top == TRUE)

  ## extract target doi peaks that must be added to tmp
  doipeaks <- target %>%
    filter(peakgr == matpkg$target_peakgr)

  ## extract tmp peaks that must be merged/averaged with doi peaks
  tmppeaks <- mat %>%
    filter(peakgr ==  matpkg$peakgr,
           target_peakgr == matpkg$target_peakgr)

  ## if target peakgr was only assigned to tmp peakgr because no other better matches were available, and in reality no peaks match up exactly, do not merge these peakgrs
  if (nrow(tmppeaks) == 0 ) {

    update <- addUNGROUPED(target = target, pkg = matpkg$target_peakgr, itmp = itmp, doi = doi)
    return(update)

  }
  else {
    
    max_peakid <- max(itmp$peakid) ## store maximum peakid for new peakid generation
      
    ## if previously grouped doi peaks were a worse fit and have to be removed from peakgr
    if (length(na.omit(previous$doi_peakid)) > 0) {

      ## update previously added doi peaks
      previous_peakid <- previous %>%
        filter(!is.na(doi_peakid)) %>%
        pull(doi_peakid)
      peakid_new <- (max_peakid + 1):(max_peakid + length(previous_peakid))
      peakgr_new <- max(itmp$peakgr) + 1
      previouspeaks <- doi %>% ## use original doi values, instead of averaged itmp values
        filter(peakid %in% previous_peakid) %>%
        mutate(doi_peakid = peakid, doi_peakgr = peakgr, doi_peakgrcls = peakgrcls) %>%
        mutate(peakid = peakid_new, peakgr = peakgr_new, peakgrcls = NA, chemid = NA, dbid = NA, dbname = NA, cos = NA) %>%
        select(peakid, mz, rt, into, peakgr, peakgrcls, doi_peakid, doi_peakgr, doi_peakgrcls, chemid, dbid, dbname, cos)
      
      max_peakid <- max(previouspeaks$peakid)
      
      ## update doi assignmnent info for the previous peaks
      doi[which(doi$peakid %in% previouspeaks$doi_peakid),
          c("tmp_peakid", "tmp_peakgr", "chemid", "dbid", "dbname", "cos")] <- previouspeaks %>%
        select(peakid, peakgr, chemid, dbid, dbname, cos)
    }
    else {
      ## if no previous tmp peaks have to be updated
      previouspeaks <- data.frame(stringsAsFactors = F)
    }

    ## merge current doi peaks with the selected tmp peakgr
    ## for every doi peak:
    ## 1. find the closest peak in the template (if multiple matches) and average mz and rt across tmp and doi
    ## 2. add chemical db metadata
    newpeaks <- doipeaks %>%
      group_by(peakid) %>%
      do(averagePEAKS(p = ., tmppeaks = tmppeaks, matpkg = matpkg)) %>%
      # ## if tmp peakid was matched by more than one doi peakid, merge doi peakids into one entry 
      # group_by(peakid) %>% 
      # do(mergePEAKS(p = .)) %>% 
      ungroup() %>%
      mutate(chemid = unique(tmppeaks$chemid), dbid = unique(tmppeaks$dbid), dbname = unique(tmppeaks$dbname))

    if (any(is.na(newpeaks$peakid))) {
      peakid_new <- (max_peakid + 1):(max_peakid + nrow(newpeaks[which(is.na(newpeaks$peakid)),]))
      newpeaks[which(is.na(newpeaks$peakid)),"peakid"] <- peakid_new
    }
    
    ## update doi peaks with assigned peakgr and peakid
    doi[which(doi$peakid %in% na.omit(newpeaks$doi_peakid)),
        c("tmp_peakid", "tmp_peakgr", "chemid", "dbid", "dbname", "cos")] <- newpeaks %>%
      filter(!is.na(doi_peakid)) %>%
      select(peakid, peakgr, chemid, dbid, dbname, cos)
    doi[which(doi$peakid %in% na.omit(newpeaks$doi_peakid)),"added"] <- rep(TRUE, length(na.omit(newpeaks$doi_peakid)))

    ## add new peaks and updated previous peaks (if any) into the structure of the tmp
    utmp <- bind_rows(newpeaks, previouspeaks)
    utmp <- bind_rows(utmp,
                      previous %>%
                        filter(!peakid %in% newpeaks$peakid) %>% ## add previous tmp peaks that were originally in the tmp and did not have to be updated
                        filter(is.na(doi_peakid)) %>%
                        mutate(doi_peakgr = unique(newpeaks$doi_peakgr),
                               doi_peakgrcls = unique(newpeaks$doi_peakgrcls),
                               cos = unique(newpeaks$cos)))

    return(list("doi" = doi, "utmp" = utmp))
  }
}

###---- if component was not assigned/grouped with any template's CID, add its peaks to template as new
addUNGROUPED <- function(target, pkg, itmp, doi ) {

  ## extract target doi peaks that must be added to tmp
  doipeaks <- target %>%
    filter(peakgr == pkg)

  peakid_new <- (max(itmp$peakid) + 1):(max(itmp$peakid) + nrow(doipeaks))
  peakgr_new <- max(itmp$peakgr) + 1
  utmp <- doipeaks %>%
    mutate(doi_peakid = peakid, doi_peakgr = peakgr, doi_peakgrcls = peakgrcls, cos = NA) %>%
    mutate(peakid = peakid_new, peakgr = peakgr_new, peakgrcls = NA) %>%
    mutate(chemid = NA, dbid = NA, dbname = NA) %>%
    select(peakid, mz, rt, into, peakgr, chemid, dbid, dbname, peakgrcls, doi_peakid, doi_peakgr, doi_peakgrcls, cos)

  doi[which(doi$peakid %in% utmp$doi_peakid),
      c("tmp_peakid", "tmp_peakgr", "chemid", "dbid", "dbname", "cos")] <- utmp %>%
    select(peakid, peakgr, chemid, dbid, dbname, cos)
  doi[which(doi$peakid %in% utmp$doi_peakid),"added"] <- rep(TRUE, length(utmp$doi_peakid))

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
    ungroup() %>%
    filter(dif == min(dif)) %>%
    ## if more than 1 peakid has the same MZ difference to target, take the most intense peakid
    filter(into == max(into))

  if (nrow(closest) > 0) {

    if (nrow(closest) > 1) {
      stop("averagePEAKS fails to select a single peak from multiple keys")
    } else {
      ## change tmp peak's mz,rt to average, into to the doi
      peaks <- closest %>%
        mutate(mz = median(c(mz, p$mz)), rt = median(c(rt, p$rt)), into = p$into) %>%
        mutate(doi_peakid = p$peakid, doi_peakgr = p$peakgr, doi_peakgrcls = p$peakgrcls, cos = matpkg$cos)
      }
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

mergePEAKS <- function(p) {
  
  if (nrow(p) > 1) {
    peak <- p %>% 
      mutate(doi_peakid = paste0(doi_peakid, collapse = ",")) %>% 
      filter(into == max(into))
  }
  else {
    peak <- p
  }
 
  return(peak)
  
}

