# show ------------------------------------------------------------------------------------------------------
#' @include classes.R
#' 
#' @rdname massFlowDB-class
#' 
#' @export
#'
setMethod("show", signature = "massFlowDB", function(object) {
  cat("A \"massFlowDB\" object \n")
  cat("Database filepath:", object@filepath, "\n")
  cat("Contains:", length(unique(object@db$dbid)), "chemicals \n")
})

# filepath --------------------------------------------------------------------------------------------------------
#' @include classes.R
#' 
#' @rdname massFlowDB-class
#' 
#' @title Obtain absolute path to chemical reference database file.
#' 
#' @description Function returns the absolute path to the \code{csv} file used when building the \code{massFlowDB} object.
#' File contains chemical reference database.
#'  
#' @export
#'
setMethod("filepath", signature = "massFlowDB", function(object) object@filepath)

# chemicals -------------------------------------------------------------------------------------------------------
#' @include classes.R
#' 
#' @rdname massFlowDB-class
#' 
#' @title Obtain names of the chemicals in the database.
#' 
#' @description Function returns names of the chemicals in the database.
#' 
#' @param object \code{massFlowDB} class object
#' 
#' @export
#'
setMethod("chemicals", signature = "massFlowDB", function(object) {
  cat("Database contains:",  unique(object@db$dbname))
})
  
