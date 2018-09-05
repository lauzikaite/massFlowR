#' Group components across peak tables.
#'
#' @param files A \code{character} with paths to peak tables with components and clusters. Created by \emph{buildCOMPS} function.
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins A \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.1).
#'
#' @return Function assigns components IDs for every component in every inputed peak table, i.e. \emph{file}. Writes updated peaks tables in the same directory as the original ones. s
#' @export
#'
#'
groupCOMPS <- function(files, mz_err = 0.01, rt_err = 2, bins = 0.1) {

  ####---- notes on default input params:
  ## bins smaller than 0.1 tend to break up really similar mzs into more bins, actually decreasing cosine accuracy
  ## rt_err of smaller than 1 tend to miss small, random variance even in almost identical tables

  ####---- create first template with mz, rt and component id for matching using the first datafile in the list
  ## save first file's table in the final output format for next stages in the pipeline:
  ## original peak table + columns:
  ## 'pid' (unique peak ID, to retain thourough grouping),
  ## 'cid' (unique component ID, to retain thourough grouping),
  ## 'clid' (cluster ID)

  message("Building first template from file: ", basename(files[[1]]))

  tmpo <- read.table(files[[1]], header = T, stringsAsFactors = F) %>%
    mutate(pid = pno, cid = comp, clid = NA, cos = NA) %>%
    select(pno, mz, rt, into, scpos, comp, pid, cid, clid, cos)
  write.table(tmpo, file = gsub(".txt", "-cid.txt" , files[[1]]), quote = F, sep = "\t", row.names = F)

  ####---- loop over all remaining datafiles in the list
  for (d in 2:length(files)) {

    message("Grouping template with file: ", basename(files[[d]]))

    ####---- update template before each new round of grouping:
    ## update tmp: remove duplicating pids, that appeared in previous grouping because of multiple doi peak matching
    ## generate new clusters for CIDs (needs to be repeated before every new doi to include newly added peaks/cids)
    ## order peaks by component complexity and peak intensity
    tmpo <- tmpo %>%
      group_by(pid) %>%
      slice(1) %>%
      ungroup()

    tmpo <- getCLUSTS(dt = tmpo %>% select(pid, mz, rt, into, cid))
    tmpo <- tmpo %>%
      select(pid, mz, rt, into, cid, clid) %>%
      mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err,
             tmp = NA, comp = NA, cls = NA, cos = NA, pno = NA, pid_to_remove = NA)

    ## a copy of tmp-original will be used to write all changes to it. Original version retains for matching comparisons
    tmp <- tmpo

    ####---- load datafile-of-interest
    ## DOI components are COMP, CIDs are assigned for each DOI peak
    doi_full <- read.table(files[[d]], header = T, stringsAsFactors = F) # save full version, as this one is going to be modified at the of the grouping round

    ## get clusters for DOI components
    doi <- doi_full %>%
      select(pid = pno, mz, rt, into, cid = comp)
    doi <- getCLUSTS(dt = doi)

    ## prepare DOI for grouping
    doi <- doi %>%
      select(pno = pid, mz, rt, into, comp = cid, cls = clid) %>%
      mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err)

    ####---- loop over all PEAKS in the DOI

    ## create a copy of DOI peaks, while() writes CIDs for each peak
    doi_peaks <- doi %>% select(pno) %>% mutate(cid = NA)

    while (any(is.na(doi_peaks$cid))) {

      ## take the first peak, which has not been assigned a CID yet
      p <- doi_peaks %>% filter(is.na(cid)) %>% slice(1) %>% pull(pno)

      # message("checking p:", p)

      target <- doi %>% filter(pno == p)

      ## extract all target peaks within CLUSTER
      target <- doi %>%
        filter(cls == target$cls)

      ####---- find MATCHES by mz/rt window and matches component (if matches in DOI are assigned to component(s), all of their features are assumed matching to this PEAK)

      ## matching by mz/er window
      mat <- target %>%
        group_by(comp, pno) %>%
        do(getMATCH(t = ., tmpo = tmpo))  %>%
        ungroup()

      ####---- COMPARE matches according to scenario (if no matches, will just update TMP)
      scen <- getSCEN(mat = mat)
      scen_out <- runSCEN(mat = mat, target = target, scen = scen, tmpo = tmpo, tmp = tmp, doi = doi, doi_peaks = doi_peaks, bins = bins, mz_err = mz_err, rt_err = rt_err)

      # if (any(!scen_out$doi_peaks %>%  filter(!is.na(cid)) %>% pull(cid) %in% tmpo$cid)) stop("new CID created, why?")

      # if (!is.na(scen_out$doi_peaks %>%  filter(pno == 394) %>% pull(cid))) stop("comp 171 was assigned")

      # if (nrow(scen_out$tmp %>%  filter(pno %in% target$pno)) != nrow(target)) stop("Not all target PNOs were added to tmp!")
      # if (nrow(scen_out$doi_peaks %>% filter(!is.na(cid)) %>% filter(!pno %in% scen_out$tmp$pno)) > 0) stop("PNO: ", p, " was overwritten in tmp!")

      tmp <- scen_out$tmp
      doi_peaks <- scen_out$doi_peaks


    }

    ###--- write grouping output
    if (all(which(doi_full$pno %in% tmp$pno)) == F) { stop("not all PNOs were assigned to TMP!")}

    ## output for doi is a table with columns:
    ## pno (original), mz (averaged), rt (averaged), into (original), scpos (original), comp (original), cls (original), pid (tmp), cid (tmp), clid (tmp), cos
    doi_full <- right_join(doi_full,
                           tmp %>%
                             filter(!is.na(pno)) %>%
                             select(pno, mz, rt, cos, pid, cid, clid), by = c("pno")) %>%
      mutate(mz = mz.y, rt = rt.y) %>%
      select(pno, mz, rt, into, scpos, comp, pid, cid, clid, cos)

    write.table(doi_full, file = gsub(".txt", "-cid.txt", files[[d]]), quote = F, sep = "\t", row.names = F)


    ## write intermediate template if grouping has to be re-initiated mid-way through
    # write.table(tmp, file = gsub(".txt", "-itmp.txt", files[[d]]), quote = F, sep = "\t", row.names = F)


    ## Update template to grouping output
    tmpo <- tmp

  }
  message("All files were succesfully grouped.")
  write.table(tmp, file = gsub(".txt", "-final-tmp.txt", files[[d]]), quote = F, sep = "\t", row.names = F)
  return(list("new_tmp" = tmp, "doi_peaks" = doi_peaks))

}
