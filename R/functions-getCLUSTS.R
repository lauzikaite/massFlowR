#' Obtain clusters of components for each datafile
#'
#' @param files
#' @param out_dir
#' @param write
#' @param workers
#'
#' @return
#' @export
#'
#' @examples
getCLUSTS <- function(files, out_dir, write = T, workers = 2) {

  BPPARAM = BiocParallel::MulticoreParam(workers = workers, log = F, progressbar = T)

  BiocParallel::bplapply(X = files,
                         FUN = buildCLUSTS,
                         out_dir = out_dir,
                         write = write,
                         BPPARAM = BPPARAM)

  message("Clusters of components were succesfully generated")

  }

buildCLUSTS <- function(f, out_dir, write){

  if (write == TRUE) { if (missing(out_dir)) { stop("'out_dir' must be specified!") } }

  fname <- gsub(".txt", "", basename(f))

  comp <- read.table(f, header = T, sep = "\t", stringsAsFactors = F)
  comp <- comp %>%
    arrange(comp)

  ## build cluster id matching table
  cls_ids <- data.frame(comp = unique(comp$comp), cls = NA)

  for(cid in unique(comp$comp)) {

    ## if component was not grouped yet into a cluster
    if (is.na(cls_ids[cid,"cls"])) {

      ## assign cluster ID
      cls_id <- ifelse(all(is.na(cls_ids$cls)), 1, max(na.omit(cls_ids$cls)) + 1)

      ## extract component's peaks
      cmp <- comp %>%
        filter(comp == cid)

      ## find all other peaks in the same RT region
      mat <- comp %>%
        # only search between components that were not clustered yet
        filter(comp %in% (cls_ids %>% filter(is.na(cls)) %>% select(comp) %>% pull) ) %>%
        filter(between(rt, min(cmp$rt), max(cmp$rt)))

      ## add other peaks in extracted components
      mat <- bind_rows(mat,
                       comp %>%
                         filter(comp %in% (mat %>% distinct(comp) %>% pull(comp))) %>%
                         filter(!pid  %in% (mat %>% distinct(pid) %>% pull(pid)))
                       )

      ## update cluster ID assignment table
      cls_ids[unique(mat$comp),"cls"] <- rep(cls_id, length(unique(mat$comp)))

    } else { next }

  }

  comp <- full_join(comp, cls_ids, by = c("comp")) %>%
    arrange(pid)

  if (write == T) {
    message(fname, " . Writing clusters table to txt file ...")
    write.table(comp, file = paste0(out_dir, fname, "-cls.txt"), col.names = T, quote = F, sep = "\t", row.names = F)
  }

  return(comp)
}

















clustCOMP <- function(tab){

  tab <- tab %>%
    arrange(comp)

  ## build cluster id matching table
  cls_ids <- data.frame(comp =  unique(tab$comp), cls = NA)

  ## set progress bar
  pb <- txtProgressBar(min = 0, max = max(tab$comp), style = 3)

  for(cid in unique(tab$comp)) {

    ## if component was not grouped yet into a cluster
    if(is.na(cls_ids[cid,"cls"])) {

      ## assign cluster ID
      cls_id <- ifelse(all(is.na(cls_ids$cls)), 1, max(na.omit(cls_ids$cls)) + 1)

      ## extract component's peaks
      cmp <- tab %>%
        filter(comp == cid)

      ## find all other peaks in the same RT region
      mat <- tab %>%
        # only search between components that were not clustered yet
        filter(comp %in% (cls_ids %>% filter(is.na(cls)) %>% select(comp) %>% pull) ) %>%
        # filter(between(rt, min(cmp$rt_l), max(cmp$rt_h))) %>%  # this is too inclusive
        # filter(rt %in% unique(cmp$rt)) %>% ## this may not be accurate enough since does not include in-between values
        filter(between(rt, min(cmp$rt), max(cmp$rt)))

      ## add other peaks in extracted components
      mat <- bind_rows(mat,
                       tab %>%
                         filter(comp %in% (mat %>% distinct(comp) %>% pull(comp))) %>%
                         filter(!order_ID  %in% (mat %>% distinct(order_ID) %>% pull(order_ID)))
                       )

      ## update cluster ID assignment table
      cls_ids[unique(mat$comp),"cls"] <- rep(cls_id, length(unique(mat$comp)))

    } else {next}

    ## update progress bar
    setTxtProgressBar(pb, cid)
  }

  tab <- full_join(tab,
                   cls_ids,
                   by = c("comp")) %>%
    arrange(order_ID)

  return(tab)
}

