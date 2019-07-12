#' @include classes.R
#'
#' @rdname massFlowAnno-class
#' 
#' @title Validate massFlowTemplate class object
#' 
validmassFlowAnno <- function(object) {
  msg <- character()
  ####---- basic validity
  if (class(object) != "massFlowAnno") {
    msg <- c(msg, "Object must be a 'massFlowAnno' class object")
  }
  ####---- validity basic slots used in initial class built
  if (!file.exists(object@filepath)) {
    msg <- c(msg, "Incorrect filepath for 'file' provided")
  }
  if (nrow(object@samples) > 0) {
    req_cnames <- c("filename",
                    "run_order")
    if (any(!req_cnames %in% names(object@samples))) {
      msg <- c(msg, paste0("'files' table must contain columns: ", paste0(req_cnames, collapse = ", ")))
    }
  }
  if (nrow(object@data) > 0) {
    req_cnames <- c("peakid",
                    "mz",
                    "rt",
                    "into",
                    "pcs")
    if (any(!req_cnames %in% names(object@data))) {
      msg <- c(msg, paste0("Slot 'data' was not correctly initiated or modified"))
    } 
  }
  if (nrow(object@ds) == 0) {
    msg <- c(msg, paste0("Slot 'ds' was not correctly initiated or modified"))
  } else {
    if (nrow(object@ds) != nrow(object@data)) {
      msg <- c(msg, paste0("Slot 'ds' must have the same number of rows as slot 'data'"))
    }
  }
  ####---- if database was loaded already (i.e. annotation was performed already)
  if (nrow(object@db) > 0) {
    ## (1) correct database table
    req_cnames <- c("peakid", "mz", "rt", "into", "chemid", "dbid", "dbname")
    if (any(!req_cnames %in% names(object@db))) {
      msg <- c(msg, paste0("Slot 'db' was not correctly initiated or modified"))
    }
    ## (2) correct parameters
    if (length(object@params) == 0) {
      msg <- c(msg, paste0("Slot 'params' doesn't contain parameters"))
    }
    ## (3) correct annotation table
    if (nrow(object@anno) == 0) {
      msg <- c(msg, paste0("Slot 'anno' was not correctly initiated or modified"))
    }
    ## (4) correct annotation results matrix
    if (nrow(object@mat) == 0) {
      msg <- c(msg, paste0("Slot 'mat' was not correctly initiated or modified"))
    }
  }
  if (length(msg)) {
    return(msg)
  } else {
    return(TRUE)
  }
}
