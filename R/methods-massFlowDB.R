setMethod("chemicals", signature = "massFlowDB", function(object) unique(object@db$dbname))
setMethod("file", signature = "massFlowDB", function(object) object@filepath)
setMethod("show", signature = "massFlowDB", function(object) {
  cat("A \"massFlowDB\" object \n")
  cat("Database filepath:", object@filepath, "\n")
  cat("Contains:", length(unique(object@db$dbid)), "chemicals \n")
})
