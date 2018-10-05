#' @aliases getDB
#' 
#' @title Get chemical reference database table using rda files
#' 
#' @description Function generates a chemical reference database table from rda files in the selected directory (dir).
#'
#' @param dir \code{character} with absolute path to the directory with rda files.
#' @param peakgr_thr \code{numeric} specifying intensity threshold (100\% BPI) for peaks in a chemical spectrum. If set to more than 0, peaks below will be removed.
#' @param compound_thr \code{numeric} specifying intensity threshold (100\% BPI) for peaks in all compound's chemical spectra. If set to more than 0, peaks below will be removed.
#' @param out_dir \code{character} with absolute path to the directory where database file should be written.
#' 
#' @details Each rda file must be written for each chemical compond separately and contain a \code{chem.file} list, comprising the following four entries:
#' \itemize{
#' \item \code{header} listing experimental details.
#' \item \code{id} specifying chemical compound identity.
#' \item \code{reference}
#' \item \code{analytical} listing all detected peaks, which were built into chemical spectra (i.e. peak groups).
#' }
#' 
#' @return Function returns generated database table and writes it to a csv file in the selected directory (out_dir).
#' 
#' @seealso \code{\link{massFlowDB}} class.
#' 
#' @export
#'
getDB <- function(dir, peakgr_thr = 0, compound_thr = 0, out_dir) {

  if (missing(dir)) stop("dir is required")
  if (!dir.exists(dir)) stop("provided path to dir is incorrect!")
  if (missing(out_dir)) stop("out_dir is required")
  if (!dir.exists(out_dir)) stop("provided path to out_dir is incorrect!")
  
  if (peakgr_thr != 0) message("peakgr_thr value above 0 is selected. Recommended value is 0")
  if (compound_thr != 0) message("compound_thr above 0 is selected. Recommended value is 0")

  db_files <- list.files(path = dir, pattern = "*.rda", full.names = T) # list DB compound rda files
  pb <- progress_estimated(n = length(db_files))
  # pb <- txtProgressBar(min = 0, max = length(db_files), style = 3)  # get progress bar
  
  ####---- for each DB compound
  db <- NULL
  for (chem in 1:length(db_files)) {
    
    load(file = db_files[chem])
    if (!exists("chem.file")) { stop("rda file is empty: ", db_files[chem]) }
    
    peakgrs <- length(chem.file$analytical) # determine the number of peakgroups in the compound's full spectrum
   
    ####---- for each peak group in the compound's file
    compound <- NULL
    for (gr in 1:peakgrs) {
      ## normalise all features IT within the single peak group to %BPI and remove features with IT < peakgr_thr
      pkg <- chem.file$analytical[[gr]]
      pkg_spec <- pkg$spec
      into_norm <- ((pkg_spec[,"IT"])/max(pkg_spec[,"IT"]))*100
      peakgr <- setNames(pkg_spec[which(into_norm > peakgr_thr),c("MZ", "IT")], nm = c("mz", "into"))
      
      if (nrow(peakgr) > 1) {
        peakgr$rt <- rep(pkg$rt, nrow(peakgr)) # add RT value for the peakgroup
        peakgr$peakgr <- rep(gr, nrow(peakgr)) # add peakgroup number
        compound <- rbind(compound, peakgr) # append to final compound matrix
      }
      else {
        next
      }
    }
    
    ## normalise all features IT within compound (all peak groups) to %BPI and remove spectrum features with IT < compound_thr
    into_norm <- ((compound[,"into"])/max(compound[,"into"]))*100
    compound <- compound[which(into_norm > compound_thr),]

    ## add chemical DB  details
    compound$chemid <- rep(chem, nrow(compound))
    compound$dbid <- rep(chem.file$id$NPCID , nrow(compound))
    compound$dbname <- rep(chem.file$id$name , nrow(compound))
    
    db <- rbind(db, compound) # append to final db matrix 
    rm(chem.file) # remove table from the environment, to make sure next file is not corrupted
    # setTxtProgressBar(pb, chem) # update progress bar
    pb$tick()$print()
    
  }  
  
  ## get unique peakgr numbers
  peakgr_ids <- db %>%
    dplyr::group_by_(.dots = c("chemid", "peakgr")) %>% 
    dplyr::group_indices()
  db$peakgr <- peakgr_ids
  
  ## arrange peaks and get unique peak numbers
  db <- db[base::with(db, order(chemid, peakgr, -into)),]
  db$rt <- db$rt * 60
  db$peakid <- 1:nrow(db)
  db <- db[,c("peakid", "mz", "rt", "into", "peakgr", "chemid", "dbid", "dbname")]
  write.csv(db, file.path(out_dir, "DBtemplate.csv"), quote = T, row.names = F)
    
  return(db)
}
