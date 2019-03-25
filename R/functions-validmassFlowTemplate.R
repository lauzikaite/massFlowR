#' @include classes.R
#'
#' @rdname massFlowTemplate-class
#' 
#' @title Validate massFlowTemplate class object
#' 
validmassFlowTemplate <- function(object) {
  msg <- character()
  ####---- basic validity
  if (class(object) != "massFlowTemplate") {
    msg <- c(msg, "Object must be a 'massFlowTemplate' class object")
  }
  ####---- validity basic slots used in initial class built
  if (!file.exists(object@filepath)) {
    msg <- c(msg, "Incorrect filepath for 'file' provided")
  }
  if (nrow(object@samples) > 0) {
    req_cnames <- c("filename",
                    "run_order",
                    "raw_filepath",
                    "proc_filepath",
                    "aligned_filepath")
    if (any(!req_cnames %in% names(object@samples))) {
      msg <- c(msg, paste0("'files' table must contain columns: ", paste0(req_cnames, collapse = ", ")))
    } else {
    if (any(!file.exists(object@samples$proc_filepath))) {
      msg <- c(msg, paste0("Column 'proc_filepath' contain incorrect file paths"))
    }
    }
  }
  if (nrow(object@tmp) > 0) {
    req_cnames <- c("peakid",
                    "mz",
                    "rt",
                    "into",
                    "peakgr")
    if (any(!req_cnames %in% names(object@tmp))) {
      msg <- c(msg, paste0("Slot 'tmp' was not correctly initiated or modified"))
    } 
  }
  
  ## if data for grouped samples is availbale
  if (length(object@data) > 0) {
    if (length(object@data) != length(object@samples$proc_filepath)) {
      ## if alignment was just initiated
      if (length(object@data) != length(which(object@samples$aligned))) {
        msg <- c(msg, paste0("Slot 'data' doesn't contain every sample in the experiment"))
      }
    }
  }
  
  if (length(object@params) == 0) {
    msg <- c(msg, paste0("Slot 'params' doesn't contain parameters"))
  }
  ####---- if validPEAKS was already applied
  if (nrow(object@valid) > 0) {
    peaks_validated <- TRUE
    if (length(object@peaks) != nrow(object@valid)) {
      msg <- c(msg, paste0("Slot 'peaks' doesn't contain a list for every validated peak"))
    }
    if (length(object@values) != length(which(object@samples$aligned))) {
      msg <- c(msg, paste0("Slot 'values' doesn't contain a list of peak values for every sample"))
    }
  } else {
    peaks_validated <- FALSE
  }
  if (length(msg)) {
    return(msg)
  } else {
    return(TRUE)
  }
}
