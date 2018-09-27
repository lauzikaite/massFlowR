#' Build a chemical database template
#'
#' @param file A \code{character} with absolute path to the the database template file (csv).
#'
#' @description Functions prepares a \code{massFlowDB} class object, which will be used for chromatograpic peaks annotations in study samples.
#'
#' @return A \code{massFlowDB} class object.
#'
#' @examples
#' ## Build template using a sample, provided with the package
#' data_dir <- system.file("testdata", package = "massFlowR")
#' db_fname <- file.path(data_dir, "DBtemplate.csv")
#' db <- buildDB(file = db_fname)
#'
buildDB <- function(file = NULL){

  if (is.null(file)) stop("file is required")
  if (any(!file.exists(file))) stop("incorrect filepath provided")

  object <- new("massFlowDB")
  object@filepath <- file

  ## populate slot with DB compounds
  db <- read.csv(file, header = T, stringsAsFactors = F)
  required_colnames <- c("peakid", "mz", "rt", "into", "peakgr", "chemid", "dbid", "dbname")
  if (any(!required_colnames %in% colnames(db))) {
    stop("DB file is missing columns: ",
         paste0(required_colnames[which(!required_colnames %in% colnames(db))],
                collapse = ", "))
  }

  object@db <- db
  return(object)
}

#' Build a sample alignment and annotation template
#'
#' @param file A \code{character} with path to the csv file, specifying samples filenames and their acquisition order.
#' @param db \code{NULL}, or a \code{massFlowDB} class object, built with \code{massFlowDB()} function.
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins A \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.1).
#' @param db_thrs A \code{numeric} specifying spectra similarity threshold (cosine) for first template generation with the database template (default set to 0.5).
#'
#' @description Functions builds a \code{massFlowTemplate} class object, which stores study sample information.
#'
#' @return A \code{massFlowTemplate} class object.
#'
#' @export
buildTMP <- function(file = NULL, db = NULL, mz_err = 0.01, rt_err = 2, bins = 0.01, db_thrs = 0.5) {

  if (is.null(file)) stop("'file' is required")
  if (!file.exists(file)) stop("incorrect filepath for 'file' provided")

  samples <- read.csv(file, header = T, stringsAsFactors = F)
  samples[,"aligned"] <- FALSE
  if (any(!c("filepaths", "run_order") %in% names(samples))) stop("files must contain columns 'filepaths' and 'run_order'")
  if (any(!file.exists(samples$filepaths))) stop("filepaths provided in the table 'file' are not correct")

  params <- list(mz_err = mz_err, rt_err = rt_err, bins = bins, db_thrs = db_thrs)

  object <- new("massFlowTemplate")
  object@filepath <- file
  object@db_filepath <- ifelse(is.null(db), "not used", db@filepath)
  object@samples <- samples
  object@params <- params

  doi_fname <- object@samples$filepaths[object@samples$run_order == 1]

  if (is.null(db)) {
    message(paste("db is not defined. Building template using sample:" , basename(doi_fname), "..."))
    tmp <- checkFILE(file = doi_fname)
    doi <- getCLUSTS(dt = tmp) %>%
      mutate(tmp_peakid = peakid, tmp_peakgr = peakgr, new_mz = mz, new_rt = rt, chemid = NA, dbid = NA, dbname = NA, cos = NA)
    write.csv(doi, file = gsub(".csv", "_aligned.csv", doi_fname), quote = T, row.names = F) ## quote = T to preserve complex DB entries with "-"
    tmp <- tmp %>%
      select(peakid, mz, rt, into, peakgr) %>%
      mutate(chemid = NA, dbid = NA, dbname = NA)

    } else {
      if (class(db) != "massFlowDB") stop("db must be a 'massFlowDB' class object")
      message(paste("Building template using database and sample:" , basename(doi_fname), "..."))
      out <- addDOI(tmp = db@db, doi_fname = doi_fname, mz_err = mz_err, rt_err = rt_err, bins = bins, add_db = TRUE, db_thrs = db_thrs)
      tmp <- out$tmp
      doi <- out$doi
      message("Template was succesfully built.")
    }

  object@tmp <- tmp
  object@samples[object@samples$run_order == 1,"aligned"] <- TRUE
  object@data[[doi_fname]] <- doi

  return(object)
}

