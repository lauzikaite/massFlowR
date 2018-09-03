####--- Create accesor methods
setGeneric("files", function(object) standardGeneric("files"))
setMethod("files", "massFlowTemplate", function(object) object@files)

setGeneric("show", function(object) standardGeneric("show"))
setMethod("show", "massFlowTemplate", function(object) {
  cat("An \"massFlowTemplate\" object with", nrow(object@files), "samples")
})

