####--- Setter and getter functions
## massFlowDB
setGeneric("chemicals", function(object) standardGeneric("chemicals"))
setMethod("chemicals", signature = "massFlowDB", function(object) unique(object@db$dbname))

setGeneric("files", function(object) standardGeneric("files"))
setMethod("files", signature = "massFlowDB", function(object) object@db_filepath)

setMethod("show", signature = "massFlowDB", function(object) {
  cat("A \"massFlowDB\" object \n")
  cat("Database filepath:", object@db_filepath, "\n")
  cat("Contains:", length(unique(object@db$dbid)), "chemicals \n")
})

## massFlowTemplate
setMethod("files", signature = "massFlowTemplate", function(object) object@files)

setMethod("show", signature = "massFlowTemplate", function(object) {
  cat("A \"massFlowTemplate\" object with", nrow(object@files), "study samples")
})


