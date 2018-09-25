#' Build a chemical database template.
#'
#' @param file A \code{character} with path to the the database template file (csv).
#' @description Functions prepares a \code{massFlowDB} class object, which will be used for chromatograpic peaks annotations in study samples.
#' @return A \code{massFlowDB} class object.
#' @export
#' @examples
#' ## Build template using a sample, provided with the package
#' data_dir <- system.file("testdata", package = "massFlowR")
#' db_fname <- file.path(data_dir, "DBtemplate.csv")
#' db <- buildDB(file = db_fname)
#'
buildDB <- function(file = NULL){

  if (is.null(file)) stop("file is required")
  if (any(!file.exists(file))) stop("incorrect filepath provided.\n")

  ## populate slot with DB compounds
  object <- new("massFlowDB")
  object@filepath <- file
  db <- read.csv(file, header = T, stringsAsFactors = F)
  required_colnames <- c("peakid", "mz", "rt", "into", "peakgr", "chemid", "dbid", "dbname")
  if (any(!required_colnames %in% colnames(db))) stop("DB file is missing columns: ",
                                                       paste0(required_colnames[which(!required_colnames %in% colnames(db))], collapse = ", "))
  object@db <- db
  return(object)
}

#' Build a sample alignment and annotation template.
#'
#' @param file A \code{character} with path to the csv file, specifying samples filenames and their acquisition order.
#' @param db NULL, or a \code{massFlowDB} class object, built with \code{massFlowDB()} function.
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins A \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.1).
#' @param db_thrs A \code{numeric} specifying spectra similarity threshold (cosine) for first template generation with the database template (default set to 0.5).
#'
#' @description Functions build \code{massFlowTemplate} class object, which stores study sample information.
#' @return A \code{\link{massFlowTemplate}} class object.
#' @export
#'
#' @examples
buildTMP <- function(file = NULL, db = NULL, mz_err = 0.01, rt_err = 2, bins = 0.1, db_thrs = 0.5) {

  ## check and add filepaths to template object
  if (is.null(file)) stop("'file' is required")
  if (!file.exists(file)) stop("incorrect filepath for 'file' provided")

  samples <- read.csv(file, header = T, stringsAsFactors = F)
  if (any(!c("filepaths", "run_order") %in% names(samples))) stop("files must contain columns 'filepaths' and 'run_order'")
  samples[,"aligned"] <- FALSE

  object <- new("massFlowTemplate")
  object@filepath <- file
  object@samples <- samples

  doi_fname <- object@samples$filepaths[object@samples$run_order == 1]

  if (is.null(db)) {
    message("db is not defined. \n Using 1st study sample to build template ... ")
    tmp <- read.csv(file = doi_fname, header = T, stringsAsFactors = F)
    } else {
      if (class(db) != "massFlowDB") stop("db must be a 'massFlowDB' class object.")
      message("Building template using database and sample: " , doi_fname, "... ")
      out <- addDOI(tmp = db@db, doi_fname = doi_fname, mz_err = mz_err, rt_err = rt_err, bins = bins, add_db = TRUE, db_thrs = db_thrs)
      tmp <- out$tmp
      message("Template was succesfully built. ")
    }

  object@tmp <- tmp
  object@samples[object@samples$run_order == 1,"aligned"] <- TRUE
  return(object)
}

