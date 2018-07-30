#' Group components across datafiles.
#'
#' @param files
#' @param mz_err
#' @param rt_err
#' @param bins
#'
#' @return
#' @export
#'
#' @examples
#'
groupCOMPS <- function(files, mz_err = 0.01, rt_err = 0.2, bins = 0.01) {

  ####---- creates the first template with mz, rt and component id for matching using the first datafile in the list ----

  message("Building first template from file: ", basename(files[[1]]))

  ## versionA - matching regions calculated using user defined mz_err and rt_err
  tmp <- read.table(files[[1]], header = T, stringsAsFactors = F) %>%
    select(pid = pid, mz, rt, into, cid = comp, cls) %>%
    mutate(mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err) %>%
    mutate(tmp = NA)

  ## versionB - matching regions taken from xcms outputs
  # tmp <- read.table(files[[1]], header = T, stringsAsFactors = F) %>%
  #   select(pid = pid, mz, rt, into, cid = comp, cls,mz_l = mzmin, mz_h = mzmax, rt_l = rtmin, rt_h = rtmax)

  ## order peaks by their components' complexity (components with more peaks go first)
  tmp_c <- tmp %>%
    group_by(cid) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    ungroup()

  tmp <- tmp %>%
    arrange(factor(cid, levels = tmp_c$cid), cid, pid)

  ####---- loop over all remaining datafiles in the list ----

  for(d in 2:length(files)) {

    message("Grouping template with file: ", basename(files[[d]]))

    ####---- load datafile-of-interest ----
    ## DOI components are COMP, CIDs are assigned for each DOI peak
    doi_full <- read.table(files[[d]], header = T, stringsAsFactors = F)

    ## versionA - matching regions calculated using mz_err and rt_err
    doi <- doi_full %>%
      select(pid = pid, mz, rt, into, comp, cls) %>%
      mutate(cid = NA, mz_l = mz - mz_err, mz_h = mz + mz_err, rt_l = rt - rt_err, rt_h = rt + rt_err)

    ## versionB - matching regions taken from xcms outputs
    # doi <- doi_full %>%
    #   select(pid = pid, mz, rt, into, comp, cls,
    #          mz_l = mzmin, mz_h = mzmax, rt_l = rtmin, rt_h = rtmax) %>%
    #   mutate(cid = NA)

    ## order peaks by their components' complexity (components with more peaks go first)
    doi_c <- doi %>%
      group_by(comp) %>%
      summarise(n = n()) %>%
      arrange(desc(n)) %>%
      ungroup()

    doi <- doi %>%
      arrange(factor(comp, levels = doi_c$comp), comp, pid)


    ####---- loop over all PEAKS in the DOI ----

    ## creates a copy of DOI peaks, while() writes CIDs for each peak
    pids <- doi %>% filter(is.na(cid)) %>% select(pid, cid)

    while(any(is.na(pids$cid))) {

      ## takes the first peak, which has not been assigned a CID yet
      p <- pids %>% filter(is.na(cid)) %>% slice(1) %>% pull(pid)

      print(paste("Checking pid:", p))

      target <- doi %>% filter(pid == p)

      ## extract all target peaks within CLUSTER
      target <- doi %>%
        filter(cls == target$cls)

      ####---- find MATCHES by mz/rt window and matches component (if matches in DOI are assigned to component(s), all of their features are assumed matching to this PEAK) ----
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
      #
      # if (any(target$pid %in% c(143))) { stop()}

      ####---- COMPARE matches according to scenario (if no matches, will just update TMP) ----
      scen <- getSCEN(tmat = mat)
      matout <- runSCEN(tmat = mat, scen = scen, tmp = tmp, pids = pids, bins = bins, mz_err = mz_err, rt_err = rt_err)
      tmp <- matout$tmp
      pids <- matout$pids

    }

    ###--- write CIDs in the original DOI datafile and write in the directory of the inputed datafiles
    doi_full <- full_join(doi_full, pids, by = c("pid"))
    write.table(doi_full, file = gsub(".txt", "-cid.txt" ,files[[d]]), quote = F, sep = "\t", row.names = F)


  }
  message("Components were succesfully grouped.")
  return(list("new_tmp" = tmp, "pids" = pids))

}
