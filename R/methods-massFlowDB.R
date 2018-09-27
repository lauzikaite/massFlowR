setMethod("show", signature = "massFlowDB", function(object) {
  cat("A \"massFlowDB\" object \n")
  cat("Database filepath:", object@filepath, "\n")
  cat("Contains:", length(unique(object@db$dbid)), "chemicals \n")
})


#' Get the the absolute path to the chemical database file
#'
#' @param A \code{massFlowDB} class object
#'
#' @export
#'
setMethod("filepath", signature = "massFlowDB", function(object) object@filepath)


#' Get the list of chemicals in the database
#'
#' @param A \code{massFlowDB} class object
#'
#' @export
#'
setMethod("chemicals", signature = "massFlowDB", function(object) unique(object@db$dbname))

