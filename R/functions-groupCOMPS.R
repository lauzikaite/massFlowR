#' Group components across datafiles.
#'
#' @param files A \code{character} with paths to peak tables with components and clusters. Created by \emph{buildCOMPS} function.
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 0.2 (min).
#' @param bins A \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.01).
#'
#' @return
#' @export
#'
#' @examples
#'
groupCOMPS <- function(files, mz_err = 0.01, rt_err = 0.2, bins = 0.01) {

  ####---- create first template with mz, rt and component id for matching using the first datafile in the list

  message("Building first template from file: ", basename(files[[1]]))

  ## save first file's table in the final output format for next stages in the pipeline:
  ## original peak table + columns:
  ##'pid' (unique peak ID, to retain thourough grouping), 'cid' (unique component ID, to retain thourough grouping),
  ##''clid' (cluster ID)

  tmp <- read.table(files[[1]], header = T, stringsAsFactors = F) %>%
    mutate(pid = pno, cid = comp, clid = cls) %>%
    ungroup()

  write.table(tmp, file = gsub(".txt", "-cid.txt" , files[[1]]), quote = F, sep = "\t", row.names = F)

  ## define first template with following columns:
  ## 'pid', 'mz', 'rt', 'into', 'cid', 'cpid', 'clid', 'mz_l', 'mz_h', 'rt_l', 'rt_h', 'tmp'
  ## peak matching regions calculated using user-defined mz_err and rt_err regions around the central mz and rt values
  tmp <- tmp %>%
    select(pid, mz, rt, into, cid, clid) %>%  # in first template, PID equals to PNO of the file
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
    mutate(tmp = NA) # This will be used by getMATCH function to check whether peak was in the tmp, or came from DOI

  ## order peaks by their components' complexity (components with more peaks go first),
  ## and peaks intensity (more intense peaks in the component go first)
  tmp_c <- tmp %>%
    group_by(cid) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    ungroup()

  tmp <- tmp %>%
    arrange(factor(cid, levels = tmp_c$cid), cid, pid)

  ####---- loop over all remaining datafiles in the list

  for(d in 2:length(files)) {

    message("Grouping template with file: ", basename(files[[d]]))

    ####---- load datafile-of-interest
    ## DOI components are COMP, CIDs are assigned for each DOI peak
    doi_full <- read.table(files[[d]], header = T, stringsAsFactors = F) # save full version, as this one is going to be modified at the of the grouping round

    ## define peak matching regions using user-defined mz_err and rt_err ranges
    doi <- doi_full %>%
      select(pno, mz, rt, into, comp, cls) %>%
      mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err,
             pid = NA, cid = NA, clid = NA) # pid, cid, clid will get filled during grouping

    ## order peaks by their components' complexity (components with more peaks go first)
    doi_c <- doi %>%
      group_by(comp) %>%
      summarise(n = n()) %>%
      arrange(desc(n)) %>%
      ungroup()
    doi <- doi %>%
      arrange(factor(comp, levels = doi_c$comp), comp, pno)

    ####---- loop over all PEAKS in the DOI

    ## create a copy of DOI peaks, while() writes CIDs for each peak
    doi_peaks <- doi %>% filter(is.na(cid)) %>% select(pno, pid, cid, clid)

    # while(any(is.na(pids$cid))) {
    while (any(is.na(doi_peaks$cid))) {

      ## take the first peak, which has not been assigned a CID yet
      p <- doi_peaks %>% filter(is.na(cid)) %>% slice(1) %>% pull(pno)

      message("checking p:", p)

      target <- doi %>% filter(pno == p)

      ## extract all target peaks within CLUSTER
      target <- doi %>%
        filter(cls == target$cls)

      ####---- find MATCHES by mz/rt window and matches component (if matches in DOI are assigned to component(s), all of their features are assumed matching to this PEAK)
      ## (1) adding key number for each peak in the target
      target <- target %>%
        group_by(comp) %>%
        mutate(key = row_number()) %>%
        mutate(key_max = max(key)) %>%
        mutate(key = as.character(key), key_max = as.character(key_max)) %>%
        ungroup()

      ## (2) matching by mz/er window
      mat <- target %>%
        group_by(comp, key) %>%
        do(getMATCH(t = ., tmp = tmp))  %>%
        ungroup()

      ####---- COMPARE matches according to scenario (if no matches, will just update TMP)
      scen <- getSCEN(tmat = mat)
      scen_out <- runSCEN(tmat = mat, scen = scen, tmp = tmp, doi_peaks = doi_peaks, bins = bins, mz_err = mz_err, rt_err = rt_err)
      tmp <- scen_out$tmp
      doi_peaks <- scen_out$doi_peaks

    }

    ###--- write CIDs in the original DOI datafile and write in the directory of the inputed datafiles
    doi_full <- full_join(doi_full, doi_peaks, by = c("pno"))
    write.table(doi_full, file = gsub(".txt", "-cid.txt" ,files[[d]]), quote = F, sep = "\t", row.names = F)


  }
  message("Components were succesfully grouped.")
  return(list("new_tmp" = tmp, "doi_peaks" = doi_peaks))

}
