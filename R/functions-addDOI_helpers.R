####---- Check and load file for alignment
checkFILE <- function(file = NULL) {
  if (!file.exists(file)) {
    stop("incorrect filepath for: ", file)
  }
  dt <- read.csv(file, stringsAsFactors = F)
  required_colnames <- c("peakid", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "peakgr")
  if (any(!required_colnames %in% colnames(dt))) {
    stop("incorrect file: ", file, " \n ",
         "missing columns: ", paste0(required_colnames[which(!required_colnames %in% colnames(dt))], collapse = ", "))
  }
  return(dt)
}


####---- add missing columns to dataframe, fill with NA
addCOLS <- function(dt, cnames) {
  add <- cnames[!cnames %in% names(dt)]
  if (length(add) != 0) {
    dt[add] <- NA
  }
  return(dt)
}

####---- add error windows using user-defined mz/rt values
addERRS <- function(dt, mz_err, rt_err) {
  dt[,c("mz_l",
        "mz_h",
        "rt_l",
        "rt_h")] <- c(dt$"mz" - mz_err,
                      dt$"mz" + mz_err,
                      dt$"rt" - rt_err,
                      dt$"rt" + rt_err)
  return(dt)
}


####---- get peakgroup clusters based on rt values
## requires columns 'peakgr','peakid','rt'
## returns original dt with one new column 'peakgrcls'
getCLUSTS <- function(dt){

  ## data frame to store assigned peak-group-clusters' ids
  cls <- data.frame(peakgr = unique(dt$peakgr), peakgrcls = NA, stringsAsFactors = F)

  ## for every peak group
  for (pkg in unique(dt$peakgr)) {

    ## if peak-group was not assigned to a cluster yet
    if (is.na(cls[which(cls$peakgr == pkg),"peakgrcls"])) {

      ## get unique peak-group-cluster id
      id <- ifelse(all(is.na(cls$peakgrcls)), 1, max(na.omit(cls$peakgrcls)) + 1)

      ## find all other peaks in the same RT region as the peak-group
      peakgroup <- dt[which(dt$peakgr == pkg),]
      cluster <- dt[which(dt$peakgr %in% cls[which(is.na(cls$peakgrcls)),"peakgr"]),] # only search between peak groups that were not clustered yet
      cluster <- cluster[which(dplyr::between(cluster$rt, min(peakgroup$rt), max(peakgroup$rt))),]  # use RT region defined by the min and max values

      ## add other peaks in the extracted peak-group
      cluster <- dt[which(dt$peakgr %in% cluster$peakgr),]

      ## update cluster ID assignment table
      cls[which(cls$peakgr %in% unique(cluster$peakgr)),"peakgrcls"] <- rep(id, length(unique(cluster$peakgr)))

    }
    else next
  }

  dt <- dplyr::full_join(dt, cls, by = c("peakgr"))

  ## order peak-groups by their complexity (number of peaks), but preserve increasing peak-group numbering
  pkg_order <- order(table(dt$peakgr), decreasing = T)
  pkg_levels <- factor(dt$peakgr, levels = pkg_order)
  dt <- dt[order(pkg_levels),]
  return(dt)

}


####---- find matches for a peak from one dataset in another
## use peak's and dataset's mz and rt regions only

matchPEAK <- function(peak, tmp) {

  ## find matching template peaks using mz/rt windows of both peaks
  mat <- tmp[which(dplyr::between(tmp$mz_l, peak["mz_l"], peak["mz_h"]) | dplyr::between(tmp$mz_h,  peak["mz_l"], peak["mz_h"])),]
  mat <- mat[which(dplyr::between(mat$rt_l, peak["rt_l"], peak["rt_h"]) | dplyr::between(mat$rt_h,peak["rt_l"], peak["rt_h"])),]

  if (nrow(mat) > 0) {
    mat$"target_peakid" <- peak["peakid"]
    mat$"target_peakgr" <- peak["peakgr"]
    mat$"target_peakgrcls" <- peak["peakgrcls"]
    mat <- mat[,c("target_peakid", "target_peakgr", "target_peakgrcls", names(tmp))]
    return(mat)
  }
}

####---- compare all matched template's peaks and find the top matches for every target peakgroup
getTOPmatches <- function(mat, target, tmp, bins, add_db, db_thrs) {

  if (nrow(mat) > 0) {
    ## extract all peaks, that belong to the same cluster as the matched cluster(s)
    allmat <- tmp[which(tmp$peakgrcls %in% mat$peakgrcls),]
    allmat <- dplyr::full_join(mat, allmat, by = names(allmat))
    
    ## check similarity of each target peakgroup against all its matches in the template
    ## for each peakgrcls-peakgrc combination:
    ## build final comparison table with both target and template peak-groups
    target_peakgr <- setNames(unique(target[c("peakgr", "peakgrcls")]), nm = c("target_peakgr", "target_peakgrcls"))
    target_peakgr_cos <- apply(target_peakgr, 1,
                               FUN = compareSPECTRA,
                               target = target, allmat = allmat, bins = bins)
    cos <- do.call(function(...) rbind(..., make.row.names = F), target_peakgr_cos)
  
    ## add chemical db info
    chemid <- lapply(cos$peakgr, function(peakgr) {
      unique(mat[which(mat$peakgr == peakgr),"chemid"])
    })
    cos$chemid <- unlist(chemid)
    
    ## compare all matching peakgrs in their clusters, extract top matches
    mattop <- compareCLUSTERS(cos = cos, add_db = add_db, db_thrs = db_thrs)

    } else {
      mattop <- data.frame(target_peakgr = unique(target$peakgr),
                           target_peakgrcls = unique(target$peakgrcls),
                           stringsAsFactors = F) %>%
        mutate(peakgrcls = NA, peakgr = NA, cos = NA, top = FALSE)
    }

  return(mattop)

}

####---- compare spectra of the target peakgroup and all matching template peakgroups
compareSPECTRA <- function(t, target, allmat, bins) {
  
  ## extract matching template peaks for the single doi target peakgroup
  matched_peaks <- allmat[which(allmat$target_peakgr == t["target_peakgr"]),]

  if (nrow(matched_peaks) > 0) {
    
    ## extract target peaks
    target_peaks <- target[which(target$peakgr == t["target_peakgr"]),]

    ## add matches that are not matched by target peaks directly
    matched_peaks <- allmat[which(allmat$peakgr %in% unique(matched_peaks$peakgr)),] 
    
    ## build total peak table for both the target peaks and all template matches
    matched_peaks["dt"] <- "tmp"
    target_peaks["dt"] <- "target"
    all_peaks <- rbind(target_peaks[,c("peakid", "peakgr", "peakgrcls", "mz", "into", "dt")],
                 matched_peaks[,c("peakid", "peakgr", "peakgrcls", "mz", "into", "dt")],
                 make.row.names = F)
    
    ## generate mz bins using all peaks
    breaks <- data.frame(breaks = seq(from = min(all_peaks$mz), to = max(all_peaks$mz), by = bins))
    breaks$bin <- 1:nrow(breaks)
    all_peaks$bin <- findInterval(x = all_peaks$mz, breaks$breaks)
    
    ## build target vector of the target peak-group spectrum
    target_vec <- all_peaks[which(all_peaks$dt == "target"),]
    target_vec <- dplyr::full_join(breaks,
                                   target_vec[which(!duplicated(target_vec$bin)),], ## remove duplicating duplicating bins
                                   by = c("bin"))
    target_vec <- ifelse(is.na(target_vec$into), 0, target_vec$into)
    target_vec <- scaleVEC(target_vec)
                                  
    ## get the cosine of the angle between vectors of the target and the matching template peakgroups
    cos <- unique(all_peaks[which(all_peaks$dt == "tmp"),c("peakgr", "peakgrcls")])
    cos_angle <- lapply(cos$peakgr,
                        FUN = getCOS,
                        all_peaks = all_peaks, breaks = breaks, target_vec = target_vec)
    cos$cos <- unlist(cos_angle)
   
  } else {
    cos <- data.frame(peakgrcls = NA, peakgr = NA, cos = NA) 
  }
  
  cos$target_peakgr <- t["target_peakgr"] 
  cos$target_peakgrcls <- t["target_peakgrcls"] 
  return(cos)
}


####---- Cosine estimation between target vector and matched vector peakgroups
getCOS <- function(m, all_peaks, breaks, target_vec){
  
  ## build vector from matched peaks
  matched_vec <- all_peaks[which(all_peaks$dt == "tmp" & all_peaks$peakgr == m),]
  matched_vec <- dplyr::full_join(breaks,
                                   matched_vec[which(!duplicated(matched_vec$bin)),], ## remove duplicating duplicating bins from different PIDS
                                   by = c("bin"))
  matched_vec <- ifelse(is.na(matched_vec$into), 0, matched_vec$into)
  matched_vec <- scaleVEC(matched_vec)

  ## find the cosine of the angle
  cos_angle <- (sum(target_vec * matched_vec) )  / ( (sqrt(sum(target_vec * target_vec)))  * ( sqrt(sum(matched_vec * matched_vec)) ) )

  return(cos_angle)
}


####---- Scale vector to unit length
scaleVEC <- function(x) {
  x / ( sqrt(sum(x * x)) )
}

####---- Compare clusters and decide on top-pairs
compareCLUSTERS <- function(cos, add_db, db_thrs) {

  ## if aligning 1st doi in the study with the db, use similarity threshold (db_thrs)
  if (add_db == T) {
    mattop <- cos %>%
      group_by(.data$peakgr) %>%
      mutate(db = ifelse(!is.na(.data$chemid), T, F)) %>%
      filter((.data$db == T & .data$cos > db_thrs) |
               (.data$db == F & .data$cos > 0)) %>%
      ungroup()
  } else {
    mattop <- cos %>%
      filter(.data$cos > 0) %>%
      ungroup()
  }

  if (nrow(mattop) > 0 ) {

   ## convert to ranks based on cosine
   mattop <- mattop %>%
     ## if pairs with the highest cosine have identical cosines (only possible when running with simulated datasets, to be investigated)
     add_count(.data$cos) %>%
     mutate(duplicated = ifelse(
       .data$n > 1 & .data$cos == max(.data$cos), T, F))
   if (any(mattop$duplicated == "TRUE")) {
     stop("identical cosines were found!")
   }

   ## (A) extract top pairs - each peakgr in the pair are ranked as 1 to each other
   mattop <- mattop %>%
     ## for each template and target peakgr separately, order the matches based on cosine
     group_by(.data$peakgr) %>%
     mutate(template_rank = row_number(desc(.data$cos))) %>%
     group_by(.data$target_peakgr) %>%
     mutate(target_rank = row_number(desc(.data$cos))) %>%
     ## for each target peakgr, take the top match (if available)
     # group_by(.data$target_peakgr) %>%
     mutate(top_rank =
              if (any(.data$template_rank == 1 & .data$target_rank == 1)) {
                ifelse(.data$template_rank == 1 & .data$target_rank == 1, TRUE, FALSE)
              } else {
                ## NA indicates that no TOP was assigned for this peakgr
                NA
              }) %>%
     ## for each tmp peakgr, mark if top was assigned, if not, mark top = NA
     group_by(.data$peakgr) %>%
     mutate(top = if (any(na.omit(.data$top_rank) == TRUE)) {
       ifelse(!is.na(.data$top_rank), ifelse(.data$top_rank == TRUE, TRUE, FALSE), FALSE)
       } else {
         NA
       })

   ## (B) extract next best pairs, by taking the pair with the highest cosine among the remaining tmp peakgr
   while (any(is.na(mattop$top))) {
     mattop <- maximiseASSIGNMENT(dt = mattop)
   }

   mattop <- mattop %>%
     group_by(.data$target_peakgr) %>%
     filter(.data$top == TRUE) %>%
     select(names(cos), .data$top) %>%
     ungroup()

 }

 ## add target peakgroups that did not have a match, will be needed when updating doi dataframe
 mattop <- mattop %>%
   bind_rows(.,
             cos %>%
               filter(!.data$target_peakgr %in% mattop$target_peakgr) %>%
               mutate(top = FALSE))

  return(mattop)

}


maximiseASSIGNMENT <- function(dt) {

  assigned <- dt %>%
    group_by(.data$target_peakgr) %>%
    filter(!is.na(.data$top)) %>%
    ungroup()

  unassigned <- dt %>%
    # group_by(.data$peakgr) %>%
    filter(is.na(.data$top)) %>%
    filter(!.data$peakgr %in% assigned$peakgr[assigned$top == TRUE]) %>%
    filter(!.data$target_peakgr %in% assigned$target_peakgr[assigned$top == TRUE])
  
  if (nrow(unassigned) > 0) {

  assigned_target_peakgrs <- assigned$target_peakgr
  assigned_target_peakgrs <- ifelse(length(assigned_target_peakgrs) == 0, 0, assigned_target_peakgrs)
  unassigned_cos <- unassigned$cos

  newly_assigned <- unassigned %>%
    # group_by(.data$peakgr) %>%
    mutate(top =
        ## if this peakgr has the highest cosine among the remaining peakgrs
        if (any(.data$cos == max(unassigned_cos))) {
          ## if its target_peakgr hasn't been assigned already
          ifelse(.data$cos == max(unassigned_cos), TRUE, FALSE)
        } else {
          NA
        })

  out <- bind_rows(assigned, newly_assigned)

  } else {
    out <- assigned
  }
  return(out)
}




####---- add doi peaks to tmp
addPEAKS <- function(mattop, mat, target, tmp, itmp, doi) {

  ## get vector with final column names order
  cnames <- c("peakid", "mz", "rt", "into", "peakgr", "peakgrcls",
              "chemid", "dbid", "dbname", 
              "doi_peakid", "doi_peakgr", "doi_peakgrcls", "cos")
  
  ## for every target/doi peakgr
  for(target_peakgr in unique(mattop$target_peakgr)) {

    ## 1. extract target/doi peaks that must be added to tmp
    target_peaks <- target[which(target$peakgr == target_peakgr),]
    
    ## 2. extract top matched tmp peakgr
    matched <- mattop[which(mattop$target_peakgr == target_peakgr),]
    matched_top <- matched[which(matched$top == T),]
    top_matched_peakgr <- ifelse(nrow(matched_top) > 0, matched_top$peakgr, 0)
    ## extract tmp peaks that must be merged/averaged with doi peaks
    matched_peaks <- mat[which(mat$peakgr == matched_top$peakgr &
                                 mat$target_peakgr == matched_top$target_peakgr),]
    
    ## if target peakgr was only assigned to tmp peakgr because of overall spectral similarity without exact matches
    if (nrow(matched_peaks) == 0) {
      ## add target peaks as new
      update <- addUNGROUPED(target_peaks = target_peaks, itmp = itmp, doi = doi, cnames = cnames)
    } else {
      ## 3. extract tmp peaks of the top matched tmp peakgr:
      ## these include both original tmp peaks, and doi peaks previously added to the tmp
      previous_peaks <- itmp[which(itmp$peakgr == top_matched_peakgr),]
      
      ## "cos" column with NA or non-existing column would break cos value comparison
      if (nrow(previous_peaks) > 0) {
        previous_peaks$cos[is.na(previous_peaks$cos)] <- 0 
      } else {
        previous_peaks <- dplyr::bind_rows(previous_peaks, data.frame(cos = 0, stringsAsFactors = F))
      }
  
      ## if current grouping is better than previous
      if (any(matched_top$cos > unique(previous_peaks$cos))) {
        update <- addGROUPED(matched_peaks, target_peaks, previous_peaks, matched_top = matched_top, itmp = itmp, doi = doi, cnames = cnames)
      } else {
        ## if current assignment is worse than previous, add target peaks as new
        update <- addUNGROUPED(target_peaks = target_peaks, itmp = itmp, doi = doi, cnames = cnames)
      }
    }

    ####---- universal part for updating tmp and doi
    doi <- update$doi
    utmp <- update$utmp
    ttmp1 <- itmp[which(!itmp$peakid %in% utmp$peakid),] ## save tmp copy without the tmp peaks matched by mz/rt window
    ttmp2 <- itmp[which(itmp$doi_peakid %in% na.omit(utmp$doi_peakid)),]  ## save tmp copy without old doi peaks that now are updated
    ttmp <- dplyr::bind_rows(ttmp1, ttmp2)
    itmp <- dplyr::bind_rows(ttmp, utmp)
    
  }
  return(list("doi" = doi, "itmp" = itmp))
}

####---- add doi peaks to tmp: peaks were grouped with an existing tmp peak-group
addGROUPED <- function(matched_peaks, target_peaks, previous_peaks, matched_top, itmp, doi, cnames) {

  ## get current maximum peakid for new peakid generation
  max_peakid <- max(itmp$peakid) 
    
  ####---- if previously added doi peaks had a worse cos (defined before this function):
  ## 1. remove from current peak-group and tmp
  ## 2. generate new peak-group
  ## 3. add to tmp under new peak-group
  if (length(na.omit(previous_peaks$doi_peakid)) > 0) {
    
    previous_peakid <- previous_peaks[which(!is.na(previous_peaks$doi_peakid)),]$doi_peakid
    peakid_new <- (max_peakid + 1):(max_peakid + length(previous_peakid))
    peakgr_new <- max(itmp$peakgr) + 1
    
    previous_peaks <- dplyr::bind_cols(
      ## doi_peakgr, doi_peakid, doi_peakgrcls remain the same
      previous_peaks[which(!is.na(previous_peaks$doi_peakid)),c("doi_peakid", "doi_peakgr", "doi_peakgrcls")],
      ## use original doi mz/rt/into values
      doi[which(doi$peakid %in% previous_peakid), c("mz", "rt", "into")])
    previous_peaks$peakid <- peakid_new
    previous_peaks$peakgr <- peakgr_new
    previous_peaks[,c("peakgrcls", "chemid", "chemid", "dbid", "dbname", "cos")] <- NA
    previous_peaks <- previous_peaks[,cnames]
    
    ## update max peakid for new peakid generation
    max_peakid <- max(peakid_new)

    } else {
    ## if no previous tmp peaks have to be updated
    previous_peaks <- data.frame(stringsAsFactors = F)
  }

  ####---- merge current target/doi peaks with the selected tmp peakgr
  ## for every doi peak:
  ## 1. find the closest peak in the template (if multiple matches) and average mz and rt across tmp and doi
  new_peaks <- apply(target_peaks, 1,
                     FUN = averagePEAKS,
                     matched_peaks = matched_peaks, matched_top = matched_top, cnames = cnames)
  new_peaks <- do.call(function(...) rbind(..., make.row.names = F), new_peaks)
  
  ## 2. add chemical db metadata
  new_peaks[,c("chemid", "dbid", "dbname")] <- unique(matched_peaks[,c("chemid", "dbid", "dbname")])

  ## 3. generate peakid for newly added doi peaks without exact peak matches
  if (any(is.na(new_peaks$peakid))) {
    peakid_new <- (max_peakid + 1):(max_peakid + nrow(new_peaks[which(is.na(new_peaks$peakid)),]))
    new_peaks[which(is.na(new_peaks$peakid)),"peakid"] <- peakid_new
  }
  
  ####---- update tmp
  ## 1. add new peaks and modified previously grouped doi peaks (if any)
  utmp <- dplyr::bind_rows(new_peaks, previous_peaks)
  
  ## 2. add tmp peaks that were originally in the tmp and did not have to be updated
  original_peaks <- itmp[which(itmp$peakgr == unique(matched_peaks$peakgr)),]
  original_peaks <- original_peaks[which(!original_peaks$peakid %in% new_peaks$peakid &
                                           is.na(original_peaks$doi_peakid)),]
  if (nrow(original_peaks) > 0) {
    original_peaks[,"doi_peakgr"] <- unique(new_peaks$doi_peakgr)
    original_peaks[,"doi_peakgrcls"] <- unique(new_peaks$doi_peakgrcls)
    original_peaks[,"cos"] <- unique(new_peaks$cos)
    
    utmp <-  dplyr::bind_rows(utmp, original_peaks)
  }
  
  utmp <- utmp[,cnames]
  
  ####---- update doi
  doi[which(doi$peakid %in% na.omit(utmp$doi_peakid)),
      c("tmp_peakid", "tmp_peakgr", "chemid", "dbid", "dbname", "cos")] <- utmp[which(!is.na(utmp$doi_peakid)),
                                                                                     c("peakid", "peakgr", "chemid", "dbid", "dbname", "cos")]
  doi[which(doi$peakid %in% na.omit(utmp$doi_peakid)),"added"] <- TRUE
  
  return(list("doi" = doi, "utmp" = utmp))
  }


###---- if component was not assigned/grouped with any template's CID, add its peaks to template as new
addUNGROUPED <- function(target_peaks, itmp, doi, cnames) {

  peakid_new <- (max(itmp$peakid) + 1):(max(itmp$peakid) + nrow(target_peaks))
  peakgr_new <- max(itmp$peakgr) + 1
  
  ## update tmp: add new peaks under new peakgr
  utmp <- target_peaks
  utmp[,c("doi_peakid", "doi_peakgr", "doi_peakgrcls")] <- utmp[,c("peakid", "peakgr", "peakgrcls")]
  utmp$peakid <- peakid_new
  utmp$peakgr <- peakgr_new
  utmp[,c("peakgrcls", "chemid", "dbid", "dbname", "cos")] <- NA
  utmp <- utmp[,cnames]

  ## update doi
  doi[which(doi$peakid %in% utmp$doi_peakid),
      c("tmp_peakid", "tmp_peakgr", "chemid", "dbid", "dbname", "cos")] <- utmp[which(!is.na(utmp$doi_peakid)),
                                                                                c("peakid", "peakgr", "chemid", "dbid", "dbname", "cos")]
  doi[which(doi$peakid %in% utmp$doi_peakid),"added"] <- TRUE

  return(list("doi" = doi, "utmp" = utmp))

}


####---- Select closest-matching peaks for each target peak (if multiple matches)
averagePEAKS <- function(peak, matched_peaks, matched_top, cnames){

  ## select closest tmp peak
  closest <- matched_peaks[which(matched_peaks$target_peakid == peak["peakid"]),]
  
  if (nrow(closest) > 0) {
    ## find mz difference between the target peakid and the matched tmp peakid
    closest$dif <- abs(peak["mz"] - closest$mz)
    closest <- closest[which(closest$dif == min(closest$dif)),]
    ## if more than 1 peakid has the same MZ difference to target, take the most intense peakid
    closest <- closest[which(closest$into == max(closest$into)),]
    
    if (nrow(closest) > 1) {
      stop("averagePEAKS fails to select a single peak from multiple keys")
    } else {
      ## change tmp peak's m & rt to average, into to the doi
      closest[,c("mz", "rt", "into")] <- c(median(closest$mz, peak["mz"]),
                                           median(closest$rt, peak["rt"]),
                                           peak["into"])
    }
  } else {
    
    ## if target peak doesn't have a match, add it as a new peak to the template
    closest <- dplyr::bind_rows(peak) ## convert numeric to data frame
    ## peakid will be generated later on during tmp merging
    closest[,c("peakid", "peakgr", "peakgrcls")] <- c(NA, matched_top$peakgr, matched_top$peakgrcls)
  }
  closest[,c("doi_peakid", "doi_peakgr", "doi_peakgrcls", "cos")] <- c(peak["peakid"], peak["peakgr"], peak["peakgrcls"], matched_top$cos)
  closest <- closest[,cnames]
  return(closest)
}

####---- if any peaks in the final template matched by multiple doi peaks, retain the most intense value
## return rownames to be removed from the final tmp
removeDUPS <- function(peakid, tmp) {
  dup_tmp <- tmp[which(tmp$peakid == peakid),]
  rnumber_max <- rownames(dup_tmp)[which(dup_tmp$into == max(dup_tmp$into))]
  rnumber <- as.numeric(rownames(dup_tmp)[which(rownames(dup_tmp) != rnumber_max)])
  return(rnumber)
}
