# checkFILE -------------------------------------------------------------------------------------------------------
## Check and load file for alignment
checkFILE <- function(file = NULL) {
  if (!file.exists(file)) {
    stop("incorrect filepath for: ", file)
  }
  dt <- read.csv(file, stringsAsFactors = F)
  required_colnames <-
    c("peakid",
      "mz",
      "mzmin",
      "mzmax",
      "rt",
      "rtmin",
      "rtmax",
      "into",
      "peakgr")
  if (any(!required_colnames %in% colnames(dt))) {
    stop(
      "incorrect file: ",
      file,
      " \n ",
      "missing columns: ",
      paste0(required_colnames[which(!required_colnames %in% colnames(dt))], collapse = ", ")
    )
  }
  return(dt)
}


# addCOLS ---------------------------------------------------------------------------------------------------------
## add missing columns to dataframe, fill with NA
# addCOLS <- function(dt, cnames) {
#   add <- cnames[!cnames %in% names(dt)]
#   if (length(add) != 0) {
#     dt[add] <- NA
#   }
#   return(dt)
# }


# addERRS ---------------------------------------------------------------------------------------------------------
## add error windows using user-defined mz/rt values
addERRS <- function(dt, mz_err, rt_err) {
  dt[, c("mz_l",
         "mz_h",
         "rt_l",
         "rt_h")] <- c(dt$"mz"-mz_err,
                       dt$"mz"+mz_err,
                       dt$"rt"-rt_err,
                       dt$"rt"+rt_err)
  return(dt)
}

# getCLUSTS -------------------------------------------------------------------------------------------------------
## get peakgroup clusters based on rt values
## requires columns 'peakgr','peakid','rt'
## returns original dt with one new column 'peakgrcls'
getCLUSTS <- function(dt) {
  ## data frame to store assigned peak-group-clusters' ids
  cls <-
    data.frame(
      peakgr = unique(dt$peakgr),
      peakgrcls = NA,
      stringsAsFactors = F
    )
  
  ## for every peak group
  for (pkg in unique(dt$peakgr)) {
    ## if peak-group was not assigned to a cluster yet
    if (is.na(cls[which(cls$peakgr == pkg), "peakgrcls"])) {
      ## get unique peak-group-cluster id
      id <-
        ifelse(all(is.na(cls$peakgrcls)), 1, max(na.omit(cls$peakgrcls)) + 1)
      
      ## find all other peaks in the same RT region as the peak-group
      peakgroup <- dt[which(dt$peakgr == pkg), ]
      cluster <-
        dt[which(dt$peakgr %in% cls[which(is.na(cls$peakgrcls)), "peakgr"]), ] # only search between peak groups that were not clustered yet
      cluster <-
        cluster[which(dplyr::between(cluster$rt, min(peakgroup$rt), max(peakgroup$rt))), ]  # use RT region defined by the min and max values
      
      ## add other peaks in the extracted peak-group
      cluster <- dt[which(dt$peakgr %in% cluster$peakgr), ]
      
      ## update cluster ID assignment table
      cls[which(cls$peakgr %in% unique(cluster$peakgr)), "peakgrcls"] <-
        rep(id, length(unique(cluster$peakgr)))
      
    } else {
      next
    }
  }
  dt <- dplyr::full_join(dt, cls, by = c("peakgr"))
  dt <- orderPEAKS(dt)
  return(dt)
}


# orderPEAKS ------------------------------------------------------------------------------------------------------
## order peak-groups by their complexity (number of peaks), but preserve increasing peak-group numbering
orderPEAKS <- function(dt) {
  pkg_order <- order(table(dt$peakgr), decreasing = T)
  pkg_levels <- factor(dt$peakgr, levels = pkg_order)
  dt <- dt[order(pkg_levels),]
  return(dt)
}


# matchPEAK -------------------------------------------------------------------------------------------------------
matchPEAK <- function(peak, tmp) {
  ## find matching template peaks using mz/rt windows of both peaks
  mat <-
    tmp[which(dplyr::between(tmp$mz_l, peak["mz_l"], peak["mz_h"]) |
                dplyr::between(tmp$mz_h,  peak["mz_l"], peak["mz_h"])),
        ## retain only the columns that you actually need
        c("peakid", "peakgr", "rt_l", "rt_h")]
  mat <-
    mat[which(dplyr::between(mat$rt_l, peak["rt_l"], peak["rt_h"]) |
                dplyr::between(mat$rt_h, peak["rt_l"], peak["rt_h"])), ]
  
  if (nrow(mat) > 0) {
    mat$"target_peakid" <- peak["peakid"]
    mat$"target_peakgr" <- peak["peakgr"]
    return(mat)
  }
}


# getTOPmatches ---------------------------------------------------------------------------------------------------
## compare all matched template's peaks and find the top matches for every target peakgroup
getTOPmatches <- function(mat, target, tmp, bins) {
  if (nrow(mat) > 0) {
    ## compare spectra for each target_peakgr with matches (instead for each unique combination)
    target_peakgrs <- unique(mat$target_peakgr)
    cos <- lapply(
      target_peakgrs,
      FUN = getCOS,
      target = target,
      mat = mat,
      tmp = tmp,
      bins = bins
    )
    cos <- do.call("rbindCLEAN", cos)
    
    ## compare all matching peakgrs in their clusters, extract top matches
    costop <-
      compareCOS(cos = cos)
    
    ## add target peakgroups that did not have a match, will be needed when updating doi dataframe
    target_peakgrs_all <- unique(target$peakgr)
    other <-
      target_peakgrs_all[(!target_peakgrs_all %in% costop$target_peakgr)]
    if (length(other) > 0) {
      other <-
        data.frame(target_peakgr = other,
                   top = F)
      costop <- dplyr::bind_rows(costop, other)
    }
  } else {
    costop <- data.frame(
      target_peakgr = unique(target$peakgr),
      peakgr = NA,
      cos = NA,
      top = F,
      stringsAsFactors = F
    )
  }
  return(costop)
}


# getCOS ----------------------------------------------------------------------------------------------------------
getCOS <- function(t_peakgr, target, mat, tmp, bins) {
  ## extract peaks from the target peakgr
  target_peaks <-
    target[which(target$peakgr == t_peakgr), c("mz", "into")]
  
  ## extract matched peakgr for this target peakgr
  matched_peakgr <-
    unique(mat[which(mat$target_peakgr == t_peakgr), c("peakgr")])
  
  ## if target peakgr matches any peaks in the tmp by mz/rt
  if (length(matched_peakgr) > 0) {
    ## extract matched peakgr peaks
    matched_peaks <-
      tmp[which(tmp$peakgr %in% matched_peakgr), c("mz", "into", "peakgr")]
    
    ## build a mz values vector reprenting all target peaks and all matched peaks
    ## this allows to have the same mz space for all matched peakgrs for this target peakgr
    ## (all matches are compared in the same mz space)
    all_peaks <- c(target_peaks$mz, matched_peaks$mz)
    breaks <-
      data.frame(breaks = seq(
        from = min(all_peaks),
        to = max(all_peaks),
        by = bins
      ))
    breaks$bin <- 1:nrow(breaks)
    
    ## build target vector of the target peakgr spectrum
    target_vec <- buildVEC(breaks = breaks, peaks = target_peaks)
    
    ## compare target peakgr and EVERY matched peakgr
    matched_peakgrs <- unique(matched_peaks$peakgr)
    cos <- lapply(matched_peakgrs, function(m_peakgr) {
      matched_peaks <-
        matched_peaks[which(matched_peaks$peakgr == m_peakgr), c("mz", "into")]
      
      ## build vector from matched peaks
      matched_vec <-
        buildVEC(breaks = breaks, peaks = matched_peaks)
      
      ## find the cosine of the angle between peakgrs
      cos_angle <-
        (sum(target_vec * matched_vec))  / ((sqrt(sum(
          target_vec * target_vec
        )))  * (sqrt(sum(
          matched_vec * matched_vec
        ))))
      return(cos_angle)
    })
    
    cos <- data.frame(target_peakgr = t_peakgr,
                      peakgr = matched_peakgrs,
                      cos = unlist(cos))
    return(cos)
  }
}

# buildVEC --------------------------------------------------------------------------------------------------------
buildVEC <- function(breaks, peaks) {
  peaks$bin <- findInterval(x = peaks$mz, breaks$breaks)
  ## remove duplicating duplicating bins
  vec <- peaks[which(!duplicated(peaks$bin)), c("bin", "into")]
  vec <- dplyr::full_join(breaks,
                          vec,
                          by = c("bin"))
  vec <- ifelse(is.na(vec$into), 0, vec$into)
  vec <- scaleVEC(vec)
  return(vec)
}

# scaleVEC --------------------------------------------------------------------------------------------------------
## Scale vector to unit length
scaleVEC <- function(x) {
  x / (sqrt(sum(x * x)))
}

# compareCOS ------------------------------------------------------------------------------------------------------
compareCOS <- function(cos) {
  ## take only positive cosines
  top <- cos[which(cos$cos > 0), ]
  
  if (nrow(top) > 0) {
    ## if pairs with the highest cosine have identical cosines
    if (nrow(top[which(top$cos == max(top$cos)),]) > 1) {
      ## investigate every case separately before making assumption in code
      stop("identical cosines were found!")
    }
    
    ## for each template peakgr (matched peakgr), order the matches based on cosine
    matched_peakgrs <- unique(top$peakgr)
    matched_rank <- lapply(matched_peakgrs, function(peakgr) {
      pairs <- cos[which(cos["peakgr"] == peakgr), ]
      top <- order(pairs["cos"], decreasing = T)
      pairs[top, "rank"] <- 1:length(top)
      return(pairs)
    })
    
    top$top <- NA
    top_t_peakgr <- NULL
    top_m_peakgr <- NULL
    ## order remaining target peakgrs by cosine
    t_cos <-
      unlist(lapply(unique(top$target_peakgr), function(t_peakgr) {
        max(top$cos[top$target_peakgr == t_peakgr])
      }))
    t_cos_order <-
      unique(top$target_peakgr)[order(t_cos, decreasing = T)]
    t_peakgrs_levels <-
      factor(top$target_peakgr, levels = t_cos_order)
    top <- top[order(t_peakgrs_levels), ]
    
    while (any(is.na(top$top))) {
      t_peakgrs <- unique(top$target_peakgr)
      
      ## for each target peakgr, take the top match (if available)
      for (t_peakgr in t_peakgrs) {
        ## order all matched peakgr by the cosine
        t_peakgr_rank <- top[top$target_peakgr == t_peakgr,]
        m_peakgrs <-
          t_peakgr_rank[order(t_peakgr_rank$cos, decreasing = T), "peakgr"]
        
        for (m_peakgr in m_peakgrs) {
          ## only if neither of the peakgrs were assigned to any top-pair already
          if (!t_peakgr %in% top_t_peakgr &
              !m_peakgr %in% top_m_peakgr) {
            ## extract the top target peakgr for this matched peakgr that was not assigned to any top pairs yet
            m_peakgr_rank <-
              matched_rank[[which(matched_peakgrs == m_peakgr)]]
            m_peakgr_rank <- m_peakgr_rank[which(
              m_peakgr_rank$target_peakgr == t_peakgr &
                !m_peakgr_rank$target_peakgr %in% top_t_peakgr
            ),]
            m_peakgr_t_peakgr <-
              m_peakgr_rank[m_peakgr_rank$rank == min(m_peakgr_rank$rank), "target_peakgr"]
            
            ## is current target peakgr the top match for the extracted top-matched-peakgr?
            if (m_peakgr_t_peakgr == t_peakgr) {
              top$top[which(top$target_peakgr == t_peakgr &
                              top$peakgr == m_peakgr)] <-
                TRUE
              top_t_peakgr <-
                c(top_t_peakgr, t_peakgr)
              top_m_peakgr <-
                c(top_m_peakgr, m_peakgr)
            }
          } else {
            top$top[which(top$target_peakgr == t_peakgr &
                            top$peakgr == m_peakgr)] <-
              FALSE
          }
        }
      }
    }
    top <- top[top$top == TRUE, ]
  }
  rownames(top) <- NULL
  return(top)
}

# addPEAKS --------------------------------------------------------------------------------------------------------
## add doi peaks to tmp
addPEAKS <- function(mattop, mat, target, tmp, itmp, doi) {
  ## get vector with final column names order
  cnames <- c(
    "peakid",
    "mz",
    "rt",
    "into",
    "peakgr",
    "peakgrcls",
    "doi_peakid",
    "doi_peakgr",
    "doi_peakgrcls",
    "cos"
  )
  
  ## for every target/doi peakgr
  for (t_peakgr in unique(mattop$target_peakgr)) {
    ## 1. extract target/doi peaks that must be added to tmp
    target_peaks <- target[which(target$peakgr == t_peakgr), ]
    
    ## 2. extract top matched tmp peakgr
    matched <-
      mattop[which(mattop$target_peakgr == t_peakgr &
                     mattop$top == T), ]
    top_matched_peakgr <-
      ifelse(nrow(matched) > 0, matched$peakgr, 0)
    
    ## extract tmp peaks that must be merged/averaged with doi peaks
    t_peakgr_mat <- mat[which(mat$peakgr == matched$peakgr &
                                mat$target_peakgr == matched$target_peakgr), ]
    matched_peaks <-
      tmp[tmp$peakid %in% t_peakgr_mat$peakid, ]
    matched_peaks$target_peakid <- t_peakgr_mat$target_peakid
    
    ## if target peakgr was only assigned to tmp peakgr because of overall spectral similarity without exact matches
    if (nrow(matched_peaks) == 0) {
      ## add target peaks as new
      update <-
        addUNGROUPED(
          target_peaks = target_peaks,
          itmp = itmp,
          doi = doi,
          cnames = cnames
        )
    } else {
      ## 3. extract tmp peaks of the top matched tmp peakgr:
      ## these include both original tmp peaks, and doi peaks previously added to the tmp
      previous_peaks <-
        itmp[which(itmp$peakgr == top_matched_peakgr), ]
      
      ## "cos" column with NA or non-existing column would break cos value comparison
      if (nrow(previous_peaks) > 0) {
        previous_peaks$cos[is.na(previous_peaks$cos)] <- 0
      } else {
        previous_peaks <-
          dplyr::bind_rows(previous_peaks,
                           data.frame(cos = 0, stringsAsFactors = F))
      }
      
      ## if current grouping is better than previous
      if (any(matched$cos > unique(previous_peaks$cos))) {
        update <-
          addGROUPED(
            matched_peaks,
            target_peaks,
            previous_peaks,
            matched = matched,
            itmp = itmp,
            doi = doi,
            cnames = cnames
          )
      } else {
        ## if current assignment is worse than previous, add target peaks as new
        update <-
          addUNGROUPED(
            target_peaks = target_peaks,
            itmp = itmp,
            doi = doi,
            cnames = cnames
          )
      }
    }
    
    ####---- universal part for updating tmp and doi
    doi <- update$doi
    utmp <- update$utmp
    ttmp1 <-
      itmp[which(!itmp$peakid %in% utmp$peakid), ] ## save tmp copy without the tmp peaks matched by mz/rt window
    ttmp2 <-
      ttmp1[which(!ttmp1$doi_peakid %in% na.omit(utmp$doi_peakid)), ]  ## save tmp copy without old doi peaks that now are updated
    itmp <- dplyr::bind_rows(ttmp2, utmp)
    
  }
  return(list("doi" = doi, "itmp" = itmp))
}


# addGROUPED ------------------------------------------------------------------------------------------------------
## add doi peaks to tmp: peaks were grouped with an existing tmp peak-group
addGROUPED <-
  function(matched_peaks,
           target_peaks,
           previous_peaks,
           matched,
           itmp,
           doi,
           cnames) {
    ## get current maximum peakid for new peakid generation
    max_peakid <- max(itmp$peakid)
    
    ####---- if previously added doi peaks had a worse cos (defined before this function):
    ## 1. remove from current peak-group and tmp
    ## 2. generate new peak-group
    ## 3. add to tmp under new peak-group
    if (length(na.omit(previous_peaks$doi_peakid)) > 0) {
      ## extract peaks coming from previous doi peakgr: doi_peakgr, doi_peakid, doi_peakgrcls remain the same
      prev_doi_peaks <-
        previous_peaks[which(!is.na(previous_peaks$doi_peakid)), c("doi_peakid", "doi_peakgr", "doi_peakgrcls")]
      prev_doi_peakid <- prev_doi_peaks$doi_peakid
      prev_doi_peakid_new <-
        (max_peakid + 1):(max_peakid + length(prev_doi_peakid))
      prev_doi_peakgr_new <- max(itmp$peakgr) + 1
      
      prev_doi_peaks <- dplyr::bind_cols(prev_doi_peaks,
                                         ## use original doi mz/rt/into values
                                         doi[which(doi$peakid %in% prev_doi_peakid), c("mz", "rt", "into")])
      
      ## check if this correctly take-up all peaks from the peakgroup - sanity check
      prev_doi_peakgr <-
        doi[which(doi$peakgr == unique(prev_doi_peaks$doi_peakgr)), "peakid"]
      if (!all(prev_doi_peaks$doi_peakid %in% prev_doi_peakgr)) {
        stop("check prev_doi_peaks extraction in addGROUPED()")
      }
      
      prev_doi_peaks$peakid <- prev_doi_peakid_new
      prev_doi_peaks$peakgr <- prev_doi_peakgr_new
      prev_doi_peaks[, c("peakgrcls", "cos")] <-
        NA
      prev_doi_peaks <- prev_doi_peaks[, cnames]
      
      ## update max peakid for new peakid generation
      max_peakid <- max(prev_doi_peakid_new)
      
    } else {
      ## if no previous tmp peaks have to be updated
      prev_doi_peaks <- data.frame(stringsAsFactors = F)
    }
    
    ####---- merge current target/doi peaks with the selected tmp peakgr
    ## for every doi peak:
    ## 1. find the closest peak in the template (if multiple matches) and average mz and rt across tmp and doi
    new_doi_peaks <- apply(
      target_peaks,
      1,
      FUN = averagePEAKS,
      matched_peaks = matched_peaks,
      matched = matched,
      cnames = cnames
    )
    new_doi_peaks <- do.call("rbindCLEAN", new_doi_peaks)
    
    ## 2. generate peakid for newly added doi peaks without exact peak matches
    if (any(is.na(new_doi_peaks$peakid))) {
      peakid_new <-
        (max_peakid + 1):(max_peakid + nrow(new_doi_peaks[which(is.na(new_doi_peaks$peakid)), ]))
      new_doi_peaks[which(is.na(new_doi_peaks$peakid)), "peakid"] <-
        peakid_new
    }
    
    ####---- update tmp
    ## 1. add new peaks and modified previously grouped doi peaks (if any)
    utmp <- dplyr::bind_rows(new_doi_peaks, prev_doi_peaks)
    
    ## 2. add tmp peaks that were originally in the tmp
    original_peaks <-
      previous_peaks[which(
        !previous_peaks$peakid %in% new_doi_peaks$peakid &
          is.na(previous_peaks$doi_peakid)
      ), ]
    if (nrow(original_peaks) > 0) {
      original_peaks[, "doi_peakgr"] <- unique(new_doi_peaks$doi_peakgr)
      original_peaks[, "doi_peakgrcls"] <-
        unique(new_doi_peaks$doi_peakgrcls)
      original_peaks[, "cos"] <- unique(new_doi_peaks$cos)
      
      utmp <-  dplyr::bind_rows(utmp, original_peaks)
    }
    
    utmp <- utmp[, cnames]
    
    ####---- update doi
    doinames <-
      c("peakid",
        "peakgr",
        "tmp_peakid",
        "tmp_peakgr",
        "cos")
    udoi <- setNames(utmp[which(!is.na(utmp$doi_peakid)),
                          c("doi_peakid",
                            "doi_peakgr",
                            "peakid",
                            "peakgr",
                            "cos")],
                     nm = doinames)
    doi[match(udoi$peakid, doi$peakid), doinames] <-
      udoi[, doinames]
    doi[match(udoi$peakid, doi$peakid), "added"] <- TRUE
    
    return(list("doi" = doi, "utmp" = utmp))
  }


# addUNGROUPED ----------------------------------------------------------------------------------------------------
## if component was not assigned/grouped with any template's CID, add its peaks to template as new
addUNGROUPED <- function(target_peaks, itmp, doi, cnames) {
  peakid_new <-
    (max(itmp$peakid) + 1):(max(itmp$peakid) + nrow(target_peaks))
  peakgr_new <- max(itmp$peakgr) + 1
  
  ## update tmp: add new peaks under new peakgr
  utmp <- target_peaks
  utmp[, c("doi_peakid", "doi_peakgr", "doi_peakgrcls")] <-
    utmp[, c("peakid", "peakgr", "peakgrcls")]
  utmp$peakid <- peakid_new
  utmp$peakgr <- peakgr_new
  utmp[, c("peakgrcls", "cos")] <- NA
  utmp <- utmp[, cnames]
  
  ## update doi
  doi[match(na.omit(utmp$doi_peakid), doi$peakid),
      c("tmp_peakid", "tmp_peakgr", "cos")] <-
    utmp[which(!is.na(utmp$doi_peakid)),
         c("peakid", "peakgr", "cos")]
  doi[match(na.omit(utmp$doi_peakid), doi$peakid), "added"] <- TRUE
  
  return(list("doi" = doi, "utmp" = utmp))
  
}

# averagePEAKS ----------------------------------------------------------------------------------------------------
## Select closest-matching peaks for each target peak (if multiple matches)
averagePEAKS <- function(peak, matched_peaks, matched, cnames) {
  ## select closest tmp peak
  closest <-
    matched_peaks[which(matched_peaks$target_peakid == peak["peakid"]), ]
  
  if (nrow(closest) > 0) {
    ## find mz difference between the target peakid and the matched tmp peakid
    closest$dif <- abs(peak["mz"] - closest$mz)
    closest <- closest[which(closest$dif == min(closest$dif)), ]
    ## if more than 1 peakid has the same MZ difference to target, take the most intense peakid
    closest <- closest[which(closest$into == max(closest$into)), ]
    
    if (nrow(closest) > 1) {
      stop("averagePEAKS fails to select a single peak from multiple keys")
    } else {
      ## change tmp peak's m & rt to average, into to the doi
      closest[, c("mz", "rt", "into")] <-
        c(median(closest$mz, peak["mz"]),
          median(closest$rt, peak["rt"]),
          peak["into"])
    }
  } else {
    ## if target peak doesn't have a match, add it as a new peak to the template
    closest <-
      ## convert numeric(created by apply) to data frame
      dplyr::bind_rows(peak) 
    ## peakid will be generated later during tmp merging
    closest[, c("peakid", "peakgr", "peakgrcls")] <-
      c(NA, matched$peakgr, NA)
  }
  closest[, c("doi_peakid", "doi_peakgr", "doi_peakgrcls", "cos")] <-
    c(peak["peakid"], peak["peakgr"], peak["peakgrcls"], matched$cos)
  closest <- closest[, cnames]
  return(closest)
}
