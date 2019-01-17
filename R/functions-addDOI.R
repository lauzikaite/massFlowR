#' @title Align and add peaks from a sample to the template
#'
#' @description Function is used by method \code{\link{alignPEAKS}}.
#' Function aligns all peaks in a datafile-of-interest (DOI) with the template (TMP).
#' Many of the required parameters are identical to the constructor function \code{\link{buildTMP}}.
#'
#' @param tmp \code{data.frame} with latest template version.
#' @param tmp_fname \code{character} specifying absolute path where updated template will be written to.
#' @param doi_fname \code{character} specifying absolute path of the datafile-of-interest csv file.
#' @param mz_err \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.1).
#' @param write_int \code{logical} whether resulting template should be written to a separate file. Default set to FALSE, so that template is over-write to the same filename.
#'
#' @seealso \code{\link{alignPEAKS}} method.
#'
#' @return Function returns updated TMP and DOI.
#' Both TMP and DOI are written to a csv file in the same directory as the original peak table for the sample.
#'
addDOI <-
  function(tmp,
           tmp_fname,
           doi_fname,
           mz_err,
           rt_err,
           bins,
           write_int) {
    ## (A) prepare tmp for grouping with the doi where doi peaks will be added to
    ## get peakgr clusters based on their retention time, order them by complexity
    tmp <- getCLUSTS(dt = tmp)
    tmp <- addERRS(dt = tmp,
                   mz_err = mz_err,
                   rt_err = rt_err)
    tmp[, c("doi_peakid", "doi_peakgr", "doi_peakgrcls", "cos")] <-
      NA
  
    ## dataframe for storing intermediate alignment results
    itmp <- tmp
    
    ## (B) prepare doi for matching against template
    doi_full <- checkFILE(file = doi_fname)
    doi_full <- getCLUSTS(dt = doi_full)
    
    ## dataframe for storing intermediate alignment results
    doi <- addERRS(dt = doi_full,
                   mz_err = mz_err,
                   rt_err = rt_err)
    doi[, c("tmp_peakgr", "tmp_peakid", "cos")] <- NA
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
        "cos",
        "added"
      )]
    
    ## initiate progress bar
    pb <-
      dplyr::progress_estimated(n = length(unique(doi$peakgrcls)))
    
    ## (C) for every doi peak, find matching tmp peaks (if any)
    while (any(doi$added == F)) {
      ## get the first peak among the un-annotated peaks
      p <- doi[doi$added == FALSE, "peakid"][1]
      
      ## get all peaks from the same cluster as that peak
      target <- doi[doi$peakid == p, ]
      target <- doi[doi$peakgrcls == target$peakgrcls, ]
      
      ## sanity check
      if (any(!is.na(target$cos))) {
        stop("check peak-group-cluster selection!")
      }
      ## match template peaks by mz/rt window
      mat <- apply(target, 1, FUN = matchPEAK, tmp = tmp)
      mat <-
        if (is.null(mat)) {
          ## if none matches in the template were found
          data.frame()
        } else {
          do.call("rbindCLEAN", mat)
        }
      # extract top matches for every target peakgr
      mattop <-
        getTOPmatches(
          mat = mat,
          target = target,
          tmp = tmp,
          bins = bins
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
    ## add empty row to align new messages below the progress bar
    message("\n")
    
    doi_out <-
      setNames(itmp[which(!is.na(itmp$doi_peakid)), c("doi_peakid",
                                                      "peakid",
                                                      "peakgr",
                                                      "mz",
                                                      "rt",
                                                      "cos")],
               nm = c(
                 "peakid",
                 "tmp_peakid",
                 "tmp_peakgr",
                 "new_mz",
                 "new_rt",
                 "cos"
               ))
    doi_out <- dplyr::full_join(doi_full, doi_out, by = "peakid")
    write.csv(
      doi_out,
      file = gsub(".csv", "_aligned.csv", doi_fname),
      quote = T,
      row.names = F
    )
    
    ## (E) Update template
    ## remove intermediate output columns not needed for next round of alignment
    tmp <-
      itmp[, c("peakid",
               "mz",
               "rt",
               "into",
               "peakgr")]
    
    ## for any tmp peaks that were matched by multiple doi peakids, retain single, most intense doi peakid
    tmp_unique <- unique(tmp[, c("mz", "rt")])
    tmp_clean <- lapply(
      1:nrow(tmp_unique),
      FUN = cleanPEAKS,
      dt_unique = tmp_unique,
      dt = tmp
    )
    tmp <- do.call("rbindCLEAN", tmp_clean)
    
    if (write_int == T) {
      ## save updated template to a unique filename for this doi
      write.csv(
        tmp,
        file = gsub(".csv", "_tmp.csv", doi_fname),
        quote = T,
        row.names = F
      )
    } else {
      ## overwrite a single template file
      write.csv(tmp,
                file = tmp_fname,
                quote = T,
                row.names = F)
    }
    return(list("tmp" = tmp, "doi" = doi_out))
  }
