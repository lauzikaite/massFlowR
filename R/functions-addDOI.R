####---- Add DOI peaks to the DB annotations table
## add DATAFILE (doi) to the the TEMPLATE (tmp)
## template could be either the master peak table, or a DATABASE (if building template for the first time)
## returns a template
addDOI <-
  function(tmp,
           tmp_fname,
           doi_fname,
           mz_err,
           rt_err,
           bins,
           add_db = FALSE,
           db_thrs = 0) {
    ## (A) prepare tmp for grouping with the doi where doi peaks will be added to
    ## get peakgr clusters based on their retention time, order them by complexity
    tmp <- getCLUSTS(dt = tmp)
    tmp <- addERRS(dt = tmp,
                   mz_err = mz_err,
                   rt_err = rt_err)
    tmp <- addCOLS(
      dt = tmp,
      c(
        "chemid",
        "dbid",
        "dbname",
        ## if not using db to build a template, add missing columns
        "doi_peakid",
        "doi_peakgr",
        "doi_peakgrcls",
        "cos"
      )
    )
    ## dataframe for storing intermediate alignment results
    itmp <- tmp
    
    ## (B) prepare doi for matching against template
    doi_full <- checkFILE(file = doi_fname)
    doi_full <- getCLUSTS(dt = doi_full)
    
    ## dataframe for storing intermediate alignment results
    doi <- addERRS(dt = doi_full,
                   mz_err = mz_err,
                   rt_err = rt_err)
    doi <- addCOLS(dt = doi,
                   c("tmp_peakgr", "tmp_peakid", "chemid", "dbid", "dbname", "cos"))
    doi$added <- FALSE
    doi <-
      doi[, c(
        "peakid",
        "mz",
        "mz_l",
        "mz_h",
        "rt",
        "rt_l",
        "rt_h",
        "into",
        "peakgr",
        "peakgrcls",
        "tmp_peakid",
        "tmp_peakgr",
        "chemid",
        "dbid",
        "dbname",
        "cos",
        "added"
      )]
    
    ## initiate progress bar
    pb <- progress_estimated(n = length(unique(doi$peakgrcls)))
    
    ## (C) for every doi peak, find matching tmp peaks (if any)
    while (any(doi$added == F)) {
      ## get the first peak among the un-annotated peaks
      p <- doi[which(doi$added == FALSE), "peakid"][1]
      
      ## get all peaks from the same cluster as that peak
      target <- doi[which(doi$peakid == p), ]
      target <- doi[which(doi$peakgrcls == target$peakgrcls), ]
      if (any(!is.na(target$cos)))
        stop("check peak-group-cluster selection!")
      
      ## match template peaks by mz/rt window
      mat <- apply(target, 1, FUN = matchPEAK, tmp = tmp)
      mat <-
        if (is.null(mat)) {
          ## if none matches in the template were found
          data.frame()
        } else {
          do.call(function(...)
            rbind(..., make.row.names = F), mat)
        }
      
      ## extract top matches for every target peakgr
      mattop <-
        getTOPmatches(
          mat = mat,
          target = target,
          tmp = tmp,
          bins = bins,
          add_db = add_db,
          db_thrs = db_thrs
        )
      
      ## add annotated/unmatched DOI peaks to template
      update <-
        addPEAKS(
          mattop = mattop,
          mat = mat,
          target = target,
          itmp = itmp,
          tmp = tmp,
          doi = doi
        )
      
      ## sanity checks
      if (update$itmp %>% filter(!is.na(doi_peakid)) %>% group_by(doi_peakid) %>% summarise(n = n()) %>% filter(n > 1) %>% nrow() > 0)
        stop("duplicated doi peakid")
      if (update$itmp %>%  group_by(peakid) %>% summarise(n = length(unique(doi_peakgr))) %>% filter(n > 1) %>% nrow() > 0)
        stop("same peakid to multiple peakgrs")
      if (any(!doi %>% filter(!is.na(tmp_peakid)) %>% pull(peakid) %in% (update$itmp %>% filter(!is.na(doi_peakid)) %>% pull(doi_peakid))))
        stop("overwritten doi peakid")
      if (any(is.na(update$itmp$peakgr)))
        stop("missing peakgr")
      
      itmp <- update$itmp
      doi <- update$doi
      
      ## update progress bar
      pb$tick()$print()
      
    }
    
    ## (D) write grouping output
    if (!all(doi$peakid %in% itmp$doi_peakid)) {
      stop("not all PNOs were assigned to TMP!")
    }
    
    ## sanity check whether peakids match up
    peakid_match <-
      unlist(lapply(unique(doi$tmp_peakgr), function(peakgr) {
        all(doi[which(doi$tmp_peakgr == peakgr), "tmp_peakid"] %in% itmp[which(itmp$peakgr == peakgr), "peakid"])
      }))
    if (any(peakid_match == F)) {
      stop("tmp_peakid and peakid don't macth up")
    }
    message("\n") ## add empty row to align new messages below the progress bar
    
    doi_out <-
      setNames(
        itmp[which(!is.na(itmp$doi_peakid)), c("doi_peakid",
                                               "peakid",
                                               "peakgr",
                                               "mz",
                                               "rt",
                                               "chemid",
                                               "dbid",
                                               "dbname",
                                               "cos")],
        nm = c(
          "peakid",
          "tmp_peakid",
          "tmp_peakgr",
          "new_mz",
          "new_rt",
          "chemid",
          "dbid",
          "dbname",
          "cos"
        )
      )
    doi_out <- dplyr::full_join(doi_full, doi_out, by = "peakid")
    write.csv(
      doi_out,
      file = gsub(".csv", "_aligned.csv", doi_fname),
      quote = T,
      row.names = F
    ) ## quote = T to preserve complex DB entries with "-"
    
    ## (E) Update template
    ## remove intermediate output columns not needed for next round of alignment
    tmp <-
      itmp[, c("peakid",
               "mz",
               "rt",
               "into",
               "peakgr",
               "chemid",
               "dbid",
               "dbname")]
    
    ## for any tmp peaks that were matched by multiple doi peakids, retain single, most intense doi peakid
    peakid_dup <- tmp$peakid[which(duplicated(tmp$peakid))]
    if (length(peakid_dup) > 0) {
      rownumber_dup <- unlist(lapply(peakid_dup,
                                     FUN = removeDUPS,
                                     tmp =  tmp))
      tmp <- tmp[-rownumber_dup, ]
    }
    write.csv(tmp,
              file = tmp_fname,
              quote = T,
              row.names = F) ## quote = T to preserve complex DB entries with "-"
    
    return(list("tmp" = tmp, "doi" = doi_out))
  }
