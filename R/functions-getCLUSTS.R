#' Obtain clusters of components for each datafile
#'
#' @param files
#' @param out_dir
#' @param bpparam
#'
#' @return
#' @export
#'
#' @examples
getCLUSTS <- function(files, out_dir, bpparam) {

  if (missing(out_dir)) { stop("'out_dir' must be specified!") }

  ## if paral workers are not defined, use the default backend
  if (missing(bpparam)) { bpparam <- bpparam() }

  BiocParallel::bplapply(X = files,
                         FUN = buildCLUSTS,
                         out_dir = out_dir,
                         BPPARAM = bpparam)

  message("Clusters of components were succesfully generated.")

  }

#' Title
#'
#' @param f
#' @param out_dir
#'
#' @return
#' @export
#'
#' @examples
buildCLUSTS <- function(f, out_dir){

  fname <- gsub(".txt", "", basename(f))
  message("Building clusters of components for file: ", fname, " ...")

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

  write.table(comp, file = paste0(out_dir, fname, "-cls.txt"), col.names = T, quote = F, sep = "\t", row.names = F)
  return(comp)

}
