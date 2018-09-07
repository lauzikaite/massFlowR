####--- Setter and getter functions
## massFlowDB
setGeneric("chemicals", function(object) standardGeneric("chemicals"))
setMethod("chemicals", signature = "massFlowDB", function(object) unique(object@db$dbname))

setGeneric("files", function(object) standardGeneric("files"))
setMethod("files", signature = "massFlowDB", function(object) object@db_filepath)

setMethod("show", signature = "massFlowDB", function(object) {
  cat("An \"massFlowDB\" object \n")
  cat("database filepath:", object@db_filepath, "\n")
  cat("contains:", length(unique(object@db$dbid)), "chemicals \n")
})

## massFlowTemplate
setMethod("files", signature = "massFlowTemplate", function(object) object@files)

setMethod("show", signature = "massFlowTemplate", function(object) {
  cat("An \"massFlowTemplate\" object with", nrow(object@files), "study samples")
})
