## constructor functions
massFlowDB <- function(files = NULL){
  
  if (is.null(files)) stop("files is required")
  if (any(!file.exists(files))) stop("incorrect filepaths in the provided [files] list are found.\n")
  
  ## populate slot with filepaths of DB compounds
  object <- new("massFlowDB")
  object@files <- data.frame(filepaths = files, stringsAsFactors = F, row.names = 1:length(files))
  
  ## populate slot with DB compounds
  dbtables <- lapply(files, function(f) { read.csv(f, header = T, stringsAsFactors = F) })
  object@db <- do.call(rbind, dbtables)
  
  return(object)
}

massFlowTemplate <- function(db = NULL, files = NULL) {
  
  ## check and add filepaths to template object
  if (is.null(files)) stop("files is required.")
  studyfiles <- read.csv(files, stringsAsFactors = F)
  if (any(!c("filepaths", "run_order") %in% names(studyfiles))) stop("files must contain columns 'filepaths'and  'run_order'")
  
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

