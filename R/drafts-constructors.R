## constructor functions
massFlowDB <- function(file = NULL){

  if (is.null(files)) stop("file is required")
  if (any(!file.exists(file))) stop("incorrect filepath in the provided.\n")

  ## populate slot with DB compounds
  object <- new("massFlowDB")
  object@db_filepath <- file
  db <- read.csv(file, header = T, stringsAsFactors = F)
  object@db <- db
  return(object)
}

massFlowTemplate <- function(files = NULL, db = NULL) {

  ## check and add filepaths to template object
  if (is.null(files)) stop("files is required.")
  studyfiles <- read.csv(files, stringsAsFactors = F)
  if (any(!c("filepaths", "run_order") %in% names(studyfiles))) stop("files must contain columns 'filepaths'and  'run_order'")

  ## create column 'grouped' where information whether sample was already grouped will be stored
  studyfiles[,"grouped"] <- FALSE
  object <- new("massFlowTemplate")
  object@files <- studyfiles

  if (is.null(db)) {
    message("db is not defined. \n Using 1st study sample to build template ... ")

    tmp <- read.csv(file = studyfiles[which(studyfiles$run_order == 1),"filepaths"], stringsAsFactors = F)
    tmp <- getCLUSTS(dt = tmp %>% rename(pid = pno, pg = comp)) %>%
                       select(pid, mz, rt, into, pg, pgc)
  } else {
    if (class(db) != "massFlowDB") stop("db must be a 'massFlowDB' class object.")
    message("Using database and 1st study sample to build template ... ")

    tmp <- annotateDOI(db = db@db, doi = object@files$filepaths[object@files$run_order == 1])

  }

  object@tmp <- tmp
  return(object)
}

