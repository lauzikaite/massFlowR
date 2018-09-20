## constructor functions
massFlowDB <- function(file = NULL){

  if (is.null(files)) stop("file is required")
  if (any(!file.exists(file))) stop("incorrect filepath provided.\n")

  ## populate slot with DB compounds
  object <- new("massFlowDB")
  object@db_filepath <- file
  db <- read.csv(file, header = T, stringsAsFactors = F)
  required_colnames <- c("peakid", "mz", "into", "peakgr", "chemid", "dbid", "dbname")
  if (any(!required_colnames %in% colnames(db))) stop("DB file is missing columns: ",
                                                       paste0(required_colnames[which(!required_colnames %in% colnames(db))], collapse = ", "))
  object@db <- db
  return(object)
}

massFlowTemplate <- function(files = NULL, db = NULL, mz_err = 0.01, rt_err = 2, bins = 0.1) {

  ## check and add filepaths to template object
  if (is.null(files)) stop("files is required.")
  if (any(!file.exists(files))) stop("incorrect filepaths provided.\n")

  studyfiles <- read.csv(files, stringsAsFactors = F)
  if (any(!c("filepaths", "run_order") %in% names(studyfiles))) stop("files must contain columns 'filepaths' and 'run_order'")

  # ## create column 'grouped' where information whether sample was already grouped will be stored
  # studyfiles[,"grouped"] <- FALSE
  object <- new("massFlowTemplate")
  object@files <- studyfiles

  if (is.null(db)) {
    message("db is not defined. \n Using 1st study sample to build template ... ")
    tmp <- read.csv(file = studyfiles[which(studyfiles$run_order == 1),"filepaths"], stringsAsFactors = F)
    } else {
      if (class(db) != "massFlowDB") stop("db must be a 'massFlowDB' class object.")
      message("Using database and 1st study sample to build template ... ")
      tmp <- addDOI(tmp = db@db, doi = object@files$filepaths[object@files$run_order == 1], mz_err = mz_err, rt_err = rt_err, bins = bins)
      }

  object@tmp <- tmp
  return(object)
}

