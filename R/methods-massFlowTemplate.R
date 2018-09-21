setMethod("file", signature = "massFlowTemplate", function(object) object@filepath)
setMethod("show", signature = "massFlowTemplate", function(object) {
  cat("A \"massFlowTemplate\" object with", nrow(object@samples), "study samples")
})
setMethod("alignSAMPLES",
          signature = "massFlowTemplate",
          function(object) {

            if (class(object) != "massFlowTemplate") stop("db must be a 'massFlowTemplate' class object.")

            return(object)


          }

)
