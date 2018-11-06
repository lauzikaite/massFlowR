#' @title Group peaks of a single LC-MS spectrum using EIC correlation
#' 
#' @description Function builds peak-groups of co-eluting peaks that are correlated to each other in a single LC-MS spectrum.
#'
#' @param pks \code{DataFrame} containing the peak table, created by \emph{pickPEAKS} function.
#' @param eic \code{list} containing extracted ion chromatograms for each peak in the \code{pks} table.
#' @param out_dir \code{character} specifying directory where output data will be saved.
#' @param fname \code{character} specifying LC-MS filename.
#' @param thr \code{numeric} defining Pearson correlation coefficient threshold, above which peak pairs will be considered as correlated.
#' @param return \code{logical} whether to return the obtained peak-group table. Set to TRUE if working with function interactivily.
#' 
#' @return Function returns an updated peak table.
#' Column \code{"peakgr"} contains peak-group ID for each peak.
#' Peaks not correlated to any of its co-eluting peaks are removed from the table.
#' Function also writes updated peak table in the specified directory.
#' 
#' @export
#'
groupPEAKSspec <- function(pks, eic, out_dir, fname, thr = thr, return = FALSE) {

  message("Building peak-groups for file: ", fname, "...")

  ## make a table to store peak-groups during grouping
  peakgroups <- pks[,c("peakid", "mz", "scpos", "into")]
  peakgroups$peakgr <- NA
  peakids <- peakgroups$peakid

  ## select POI (peak of interest) from peak table
  for (p in peakids) {
    
    ## if POI is already assigned to a group, skip to next
    if (!is.na(peakgroups[p, "peakgr"])) next 

    ## find co-eluting peaks not asigned to a peakgronent yet
    ind <- c(peakgroups[p, "scpos"] - 1, peakgroups[p, "scpos"] + 1)
    co <- peakgroups[which(is.na(peakgroups$peakgr) &
                             dplyr::between(peakgroups$scpos, ind[1], ind[2])),]
    co_ind <- co$peakid
    
    ## if POI doesn't have co-eluting peaks, skip to next
    if (length(co_ind) == 1 | !p %in% co_ind) next 

    ## correlate the EIC of the peaks in the correlation matrix (cormat)
    cormat <- getCORmat(ind = co_ind)
    cormat$weight <- apply(cormat, 1, FUN = corEIC, eic = eic)
    
    ## if there isn't a single pair with cor above threshold, skip to next
    if (all(cormat$weight < thr)) next 
    
    ## build network between co-eluting peaks and obtain communities
    coms <- buildGRAPH(pkg_cor = cormat, cor_thr = thr, plot = FALSE)
    
    ## extract community which includes the POI
    main_com <- coms[which(names(coms) == p)]
    
    ## extract all co-eluting peaks that are part of this community
    co_ind_com <- as.numeric(names(coms[which(coms == main_com)]))
   
    ## update table with assigned peakgroup IDs
    ## start with 1 if this is the first peakgroup to be assigned
    peakid <- ifelse(all(is.na(peakgroups$peakgr)), 1, max(peakgroups$peakgr, na.rm = T) + 1)
    peakgroups[co_ind_com, "peakgr"] <- peakid
  }
  
  ## removing peaks not assigned to any peak-group
  peakgroups <- peakgroups[which(!is.na(peakgroups$peakgr)),]
  peakgroups <- merge(pks, 
                      peakgroups[, c("peakid", "peakgr")],
                    by = c("peakid"), all = F)
  ## update peakids
  peakgroups$peakid <- 1:nrow(peakgroups)
  write.csv(peakgroups,
            file = paste0(out_dir, "/", fname, "_peakgrs.csv"),
            quote = F,
            row.names = F)

  if (return == TRUE) {
    return(peakgroups)
  } else {
    message(max(na.omit(peakgroups$peakgr)), " peak-groups built.")
  }
}

# corEIC ------------------------------------------------------------------------------------------------------
#' @title Obtain EIC correlation between two peaks
#' 
#' @description Function performs extracted ion chromatogram (EIC) correlation of two peaks.
#'
#' @param pair \code{matrix} with columns 'from' and 'to', specifying peak identifiers, which correspond to the indeces of a list of EICs.
#' @param eic \code{list} containing EIC for each picked-peak in a sample. 
#'
#' @return Function returns EIC correlation coefficient.
#'
corEIC <- function(pair, eic) {
  x <- pair["from"]
  y <- pair["to"]
  rx <- MSnbase::rtime(eic[[x]])
  ry <- MSnbase::rtime(eic[[y]])
  common_scan <- base::intersect(rx, ry)
  if (length(common_scan) > 3) {
    ix <- as.numeric(MSnbase::intensity(eic[[x]])[which(rx %in% common_scan)])
    iy <- as.numeric(MSnbase::intensity(eic[[y]])[which(ry %in% common_scan)])
    cc <- cor(ix, iy, method = "pearson")
  } else {
    cc <- 0
  }
  ## negative coeficients would break graph generation
  cc <- ifelse(cc < 0, 0, cc) 
  return(cc)
}