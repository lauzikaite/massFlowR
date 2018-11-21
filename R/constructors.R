# buildDB ------------------------------------------------------------------------------------------------------
#' @aliases buildDB
#'
#' @title Build a chemical database template
#'
#' @description Function creates a \code{massFlowDB} class object, which will be used for chromatograpic peaks annotation and study samples alignment.
#'
#' @param file \code{character} with absolute path to the the database template file (csv).
#'
#' @return Function returns a \code{massFlowDB} class object.
#'
#' @export
#'
#' @seealso \code{\link{massFlowDB}} class
#'
#' @examples
#' ## Build database template using a sample file, provided with the package
#' data_dir <- system.file("testdata", package = "massFlowR")
#' db_fname <- file.path(data_dir, "DBtemplate.csv")
#' db <- buildDB(file = db_fname)
#'
buildDB <- function(file = NULL) {
  if (is.null(file)) {
    stop("file is required")
  }
  if (any(!file.exists(file))) {
    stop("incorrect filepath provided: ", file)
  }
  
  object <- new("massFlowDB")
  object@filepath <- file
  
  ## populate slot with DB compounds
  db <- read.csv(file, header = T, stringsAsFactors = F)
  required_colnames <-
    c("peakid",
      "mz",
      "rt",
      "into",
      "peakgr",
      "chemid",
      "dbid",
      "dbname")
  if (any(!required_colnames %in% colnames(db))) {
    stop("DB file is missing columns: ",
         paste0(required_colnames[which(!required_colnames %in% colnames(db))],
                collapse = ", "))
  }
  
  object@db <- db
  return(object)
}


# buildTMP --------------------------------------------------------------------------------------------------------
#' Build a sample alignment and annotation template
#'
#' Functions builds a \code{massFlowTemplate} class object, which stores study sample information.
#'
#' @param file \code{character} for absolute path to the csv file, specifying samples filenames and their acquisition order.
#' @param mz_err \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.01).
#'
#' @return A \code{massFlowTemplate} class object.
#'
#' @export
buildTMP <-
  function(file = NULL,
           mz_err = 0.01,
           rt_err = 2,
           bins = 0.01
           ) {
    if (is.null(file)) {
      stop("'file' is required")
    }
    if (!file.exists(file)) {
      stop("incorrect filepath for 'file' provided")
    }
    samples <- read.csv(file, header = T, stringsAsFactors = F)
    if (any(!c("filepaths", "run_order") %in% names(samples))) {
      stop("files must contain columns 'filepaths' and 'run_order'")
    }
    ## stop only if none of the filepaths exist - some may still being processed on another machine
    if (all(!file.exists(samples$filepaths))) {
      stop("filepaths provided in the table 'file' are not correct")
    }
    
    object <- new("massFlowTemplate")
    object@filepath <- file
    object@samples <- samples
    object@samples[, "aligned"] <- FALSE
    object@params <-
      list(mz_err = mz_err,
           rt_err = rt_err,
           bins = bins)
    tmp_fname <-
      gsub(".csv", "_template.csv", object@filepath)
    doi_first <- min(object@samples$run_order)
    doi_fname <-
      object@samples$filepaths[object@samples$run_order == doi_first]

    message(paste("Building template using sample:",
                  basename(doi_fname),
                  "..."))
    ## write 1st sample in the standard output format
    doi <- checkFILE(file = doi_fname)
    doi <- getCLUSTS(dt = doi)
    doi[, c("tmp_peakid", "tmp_peakgr", "new_mz", "new_rt")] <-
      doi[, c("peakid", "peakgr", "mz", "rt")]
    doi[,c("cos")] <- NA
    write.csv(
      doi,
      file = gsub(".csv", "_aligned.csv", doi_fname),
      ## quote = T to preserve complex DB entries with "-"
      quote = T,
      row.names = F
    ) 
    ## build template from 1st sample
    tmp <-
      doi[, c("peakid",
              "mz",
              "rt",
              "into",
              "peakgr")]

    object@tmp <- tmp
    object@samples[object@samples$run_order == doi_first, "aligned"] <-
      TRUE
    object@data[[doi_fname]] <- doi
    return(object)
  }


# loadALIGNED -----------------------------------------------------------------------------------------------------
#' Build sample alignment and annotation template using already aligned samples
#'
#' Function handles the construction of \code{massFlowTemplate} class object from already aligned samples.
#' Function enables user to continue an interrupted peak alignment process, facilitated by \code{alignPEAKS} function.
#'
#' @details Arguments are identical to the ones used by \code{\link{buildTMP}} constructor function.
#'
#' @param file A \code{character} with path to the csv file, specifying samples filenames and their acquisition order.
#' @param db \code{NULL}, or a \code{massFlowDB} class object, built with \code{massFlowDB()} function.
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins A \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.01).
#'
#' @return A \code{massFlowTemplate} class object.
#'
#' @seealso \code{\link{massFlowTemplate}} class.
#'
#' @export
#'
loadALIGNED <-
  function(file,
           db = NULL,
           mz_err = 0.01,
           rt_err = 2,
           bins = 0.01) {
    if (is.null(file)) {
      stop("'file' is required")
    }
    if (!file.exists(file)) {
      stop("incorrect filepath for 'file' provided")
    }
    samples <- read.csv(file, header = T, stringsAsFactors = F)
    if (any(!c("filepaths", "run_order") %in% names(samples))) {
      stop("files must contain columns 'filepaths' and 'run_order'")
    }
    if (all(!file.exists(samples$filepaths))) {
      stop("filepaths provided in the table 'file' are not correct")
    }
    
    ## load the latest template file
    tmp_file <- gsub(".csv", "_template.csv", file)
    if (!file.exists(tmp_file)) {
      stop("template file is not available: ", tmp_file)
    }
    filepaths_aligned <-
      gsub(".csv", "_aligned.csv", samples$filepaths)
    samples_aligned <- which(file.exists(filepaths_aligned))
    tmp <- read.csv(tmp_file, header = T, stringsAsFactors = F)
    
    object <- new("massFlowTemplate")
    object@filepath <- file
    object@samples <- samples
    object@samples[, "aligned"] <- FALSE
    object@samples[samples_aligned, "aligned"] <- TRUE
    object@params <-
      list(
        mz_err = mz_err,
        rt_err = rt_err,
        bins = bins
      )
    object@tmp <- tmp
    
    ## load aligned samples datasets
    data <-
      lapply(filepaths_aligned[samples_aligned], function(doi_fname) {
        doi <- read.csv(doi_fname,
                        header = T,
                        stringsAsFactors = F)
      })
    names(data) <- object@samples[samples_aligned, "filepaths"]
    object@data <- data
    message("A 'massFlowTemplate' object was succesfully built with aligned samples.")
    return(object)
  }