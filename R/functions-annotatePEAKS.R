#' @title Annotate dataset using a chemical reference database file
#' 
#' @description Function annotates a dataset using a chemical reference database.
#' Peak-groups, or Pseudo Chemical Spectra (PCS), in the dataset are compared with the peak-groups representing reference compounds in the database.
#' PCS are annotated using dot-product estimation between the spectra of the PCS and the corresponding chemical in the database.
#'
#' @param dataset \code{character} for absolute path to the csv file of the dataset of interest.
#' @param database \code{character} for absolute path to the csv file of the chemical reference database, built using \code{\link{buildDB}} function.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param mz_err \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 10 (sec).
#' @param bins \code{numeric} defying step size used in peak-group spectra binning and vector generation. Step size represents MZ dimension. Default set to 0.01.
#' @param ncores \code{numeric} for number of parallel workers to be used. Set 1 for serial implementation. Default set to 2.
#'
#' @return Function writes an updated dataset table with columns 'db_peakid' and 'db_chemid', indicating matched chemical reference standards from the database.
#' Column 'cos' indicates the stringth of spectral similarity between the PCS and the corresponding chemical compound, ranging from 0 to 1.
#' 
#' @export
#'
#' @seealso \code{\link{buildDB}}, \code{\link{alignPEAKS}}
#' 
annotatePEAKS <- function(dataset = NULL,
                          database = NULL,
                          out_dir = NULL,
                          mz_err = 0.1,
                          rt_err = 10,
                          bins = 0.1,
                          ncores = 2
                          ) {
  if (is.null(dataset)) {
    stop("'dataset' filepath is required")
  }
  if (is.null(database)) {
    stop("'database' filepath is required")
  }
  if (is.null(out_dir)) {
    stop("'out_dir' is required")
  }
  if (!dir.exists(out_dir)) {
    stop("incorrect filepath for 'out_dir' provided")
  }
  ## register paral backend
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  ds_original <- read.csv(dataset, header = T, stringsAsFactors = F)
  db_original <- read.csv(database, header = T, stringsAsFactors = F)
  ds_cnames <- c("peakid", "pcs", "mz", "rt", "into")
  ds <- ds_original[, c(match(ds_cnames, colnames(ds_original)))]
  db <- db_original
  ## temp fix for float precision
  ds[, c("mz", "rt", "into")] <- t(apply(ds[, c("mz", "rt", "into")], 1, round, digits = 8))
  db[, c("mz", "rt", "into")] <- t(apply(db[, c("mz", "rt", "into")], 1, round, digits = 8))
  
  ####---- compare dataset with the database
  ds_to_db <- do_alignPEAKS(
    ds = ds,
    tmp = db,
    ds_var_name = "pcs",
    tmp_var_name = "chemid",
    mz_err = mz_err,
    rt_err = rt_err,
    bins = bins,
    ncores = ncores
  )
  
  ####---- export annotation table
  ## list to store tmp peakids for every peak in doi
  ds_to_db_peakids <- rep(list(NA), nrow(ds))
  ds_to_db_chemid <- rep(list(NA), nrow(ds))
  ds_to_db_cos <- rep(list(NA), nrow(ds))
  
  ## for EVERY peak-group in dataset
  for (var in 1:length(ds_to_db)) {
    ds_var <- ds_to_db[[var]]$ds
    ds_var_peaks <- ds[ds$pcs == ds_var, ]
    db_var <-  ds_to_db[[var]]$tmp
    if (!is.null(db_var)) {
      ds_peakids <- ds_to_db[[var]]$mat$target_peakid
      db_peakids <- ds_to_db[[var]]$mat$peakid
      cos <- ds_to_db[[var]]$cos
      ## if any unmatched peaks are left
      ds_var_peaks <- ds[which(ds$pcs == ds_var), ]
      ds_var_peaks_un <- ds_var_peaks[!ds_var_peaks$peakid %in% ds_peakids, ]
      if (nrow(ds_var_peaks_un) > 0) {
        db_peakids <- c(db_peakids, rep(NA, nrow(ds_var_peaks_un)))
        ds_peakids <- c(ds_peakids, ds_var_peaks_un$peakid)
      }
    }
    ## update doi-to-tmp matching info
    ds_to_db_peakids[unlist(ds_peakids)] <- db_peakids
    ds_to_db_chemid[unlist(ds_peakids)] <- db_var
    ds_to_db_cos[unlist(ds_peakids)] <- cos
  }
  ds$db_peakid <- unlist(ds_to_db_peakids)
  ds$db_chemid <- unlist(ds_to_db_chemid)
  ds$cos <- unlist(ds_to_db_cos)
  
  write.csv(ds,
            file = file.path(out_dir, gsub(".csv", "_annotated.csv", basename(dataset))),
            quote = T,
            row.names = F)
  message("Dataset was annotated.")
}
