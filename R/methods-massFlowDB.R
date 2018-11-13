#' @aliases show
#'
#' @rdname massFlowDB-class
#' 
#' @param object \code{massFlowDB} class object
#' 
#' @export
#'
setMethod("show", signature = "massFlowDB", function(object) {
  cat("A \"massFlowDB\" object \n")
  cat("Database filepath:", object@filepath, "\n")
  cat("Contains:", length(unique(object@db$dbid)), "chemicals \n")
})


#' @aliases filepath
#' @title Get the the absolute path to the chemical database file
#'
#' @rdname massFlowDB-class
#' 
#' @param object \code{massFlowDB} class object
#' 
#' @export
#'
setMethod("filepath", signature = "massFlowDB", function(object) object@filepath)


#' Get the list of chemicals in the database
#' @aliases chemicals
#'
#' @rdname massFlowDB-class
#' 
#' @param object \code{massFlowDB} class object
#' 
#' @export
#'
setMethod("chemicals", signature = "massFlowDB", function(object) {
  cat("Database contains:",  unique(object@db$dbname))
})
  
  
 

