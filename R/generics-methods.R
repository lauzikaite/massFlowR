####--- Setter and getter functions
## massFlowDB
setGeneric("files", function(object) standardGeneric("files"))
setMethod("files", signature = "massFlowDB", function(object) object@files)

setMethod("show", "massFlowDB", function(object) {
  cat("An \"massFlowDB\" object with", nrow(object@files), "database files")
})

## massFlowTemplate
setMethod("files", signature = "massFlowTemplate", function(object) object@files)

setMethod("show", "massFlowTemplate", function(object) {
  cat("An \"massFlowTemplate\" object with", nrow(object@files), "study samples")
})