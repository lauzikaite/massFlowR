setMethod("filepath", signature = "massFlowTemplate", function(object) object@filepath)
setMethod("show", signature = "massFlowTemplate", function(object) {
  cat("A \"massFlowTemplate\" object with", nrow(object@samples), "study samples")
})
setMethod("alignSAMPLES",
          signature = "massFlowTemplate",
          function(object, mz_err = 0.01, rt_err = 2, bins = 0.01) {

            if (class(object) != "massFlowTemplate") stop("db must be a 'massFlowTemplate' class object")

            while (any(object@samples$aligned == FALSE)) {

              message(paste(length(which(object@samples$aligned == F)), "out of" , nrow(object@samples), "samples left to align."))
              doi_fname <- checkSAMPLES(object)
              message("Aligning sample: ", basename(doi_fname),  "... ")
              out <- addDOI(tmp = object@tmp, doi_fname = doi_fname, mz_err = mz_err, rt_err = rt_err, bins = bins)
              object@tmp <- out$tmp
              object@samples[object@samples$filepaths == doi_fname,"aligned"] <- TRUE
              object@data[[doi_fname]] <- out$doi

            }
            message("All samples were aligned.")
            return(object)
          }

)

setMethod("checkSAMPLES",
          signature = "massFlowTemplate",
          function(object) {
            if (class(object) != "massFlowTemplate") stop("db must be a 'massFlowTemplate' class object.")

            ## extract next-in-line sample and check if it is present (i.e. has been written already)
            doi_fname <- object@samples$filepaths[which(object@samples$aligned == FALSE)[1]]
            doi_present <- file.exists(doi_fname) # if FALSE, then loop until TRUE

            while (!doi_present) {
              if (!file.exists(doi_fname)) {
                message("Waiting for file to be written for: ", doi_fname, " ...")
                Sys.sleep(time = 2)
              }
              doi_present <- file.exists(doi_fname)
            }

            return(doi_fname)
          }
)

