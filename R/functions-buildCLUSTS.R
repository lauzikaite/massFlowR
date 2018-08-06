#' @title Component assignment into chromatographic regions (clusters).
#' @description Function assigns components of correlated peaks into chromatographic regions, called clusters.
#'
#' @param fname A \code{character} specifying datafile name.
#' @param comp_table A \code{DataFrame} class peak table with assigned components, created by \emph{buildCOMPS} function.
#' @param out_dir out_dir \code{character} object specifying directory where output data will be saved.
#'
#' @return Function returns an updated peak table. Column \code{"cls"} contains cluster ID for each peak.
#' Function also writes updated peak table to a text file in the specified directory.
#' @export
#'
#' @examples
buildCLUSTS <- function(fname, comp_table, out_dir){

  comp <- comp_table %>%
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
                         filter(!pno  %in% (mat %>% distinct(pno) %>% pull(pno)))
                       )

      ## update cluster ID assignment table
      cls_ids[unique(mat$comp),"cls"] <- rep(cls_id, length(unique(mat$comp)))

    } else { next }

  }

  cls <- full_join(comp, cls_ids, by = c("comp")) %>%
    arrange(pno)

  write.table(cls, file = paste0(out_dir, "/", fname, "_pks-comps-cls.txt"), col.names = T, quote = F, sep = "\t", row.names = F)

}
