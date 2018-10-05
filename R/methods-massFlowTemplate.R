####---- show - documentation in main massFlowTemplate-class page ----
#' @aliases show
#'
#' @rdname massFlowTemplate-class
#' 
#' @export
#' 
setMethod("show", signature = "massFlowTemplate", function(object) {
  cat("A \"massFlowTemplate\" object with", nrow(object@samples), " samples")
})

####---- filepath - documentation in main massFlowTemplate-class page ----
#' @aliases filepath
#' 
#' @rdname massFlowTemplate-class
#'
#' @export
#'
setMethod("filepath", signature = "massFlowTemplate", function(object) object@filepath)

####---- checkNEXT - no documentation, as not exported ----
setMethod("checkNEXT",
          signature = "massFlowTemplate",
          function(object) {
            if (class(object) != "massFlowTemplate") stop("db must be a 'massFlowTemplate' class object.")
            
            ## extract next-in-line sample and check if it is present (i.e. has been written already)
            doi_fname <- object@samples$filepaths[which(object@samples$aligned == FALSE)[1]]
            doi_present <- file.exists(doi_fname) # if FALSE, then loop until TRUE
            
            while (!doi_present) {
              if (!file.exists(doi_fname)) {
                message(paste("Waiting for file to be written for:", doi_fname, "..."))
                Sys.sleep(time = 60)
              }
              doi_present <- file.exists(doi_fname)
            }
            
            return(doi_fname)
          }
)


####---- alignPEAKS - more complex method, has a separate documentation page ----
#' @aliases alignPEAKS
#' 
#' @title Align peaks detected in LC-MS samples using spectral similarity comparison
#'
#' @description Method aligns peaks across samples in LC-MS experiment using spectral similarity comparison.
#' To enable alignment, peaks originating from the same chemical compound were grouped into pseudo chemical spectra (PCS), via function \code{\link{groupPEAKS}}.
#'
#' @details Peaks are aligned across samples in their original acquisition order.
#' Template is list of all previously detected and aligned peaks.
#' For each peak in a sample, \code{alignPEAKS}:
#' \itemize{
#' \item Finds all template peaks within a \emph{m/z} and \emph{rt} window.
#' \item Identifies the true match by comparing the spectral similarity between the peak-group of the peak-of-interest and all matching template's peak-groups.
#' \item Merges the selected template's peak-group with the peak-group of the peak-of-interest.
#' \item Updates template's \emph{m/z} and \emph{rt} values for the matching peaks across the template and the sample.
#' }
#'
#' Spectral similarity is measured by obtaining the cosine of the angle between two 2D vectors, representing each peak-group's \emph{m/z} and \emph{intensity} values.
#'
#' @param object \code{massFlowTemplate} class object, created by \emph{buildTMP} constructor function.
#'
#' @return Method updates \code{\link{massFlowTemplate}} class object.
#' Slot @@tmp is updated after each round of sample alignment.
#' Slot @@data stores alignment results for each sample in the experiment.
#' Method writes alignment results to a csv file for each sample in the experiment.
#' 
#' @export
#'
setMethod("alignPEAKS",
          signature = "massFlowTemplate",
          function(object) {

            if (class(object) != "massFlowTemplate") stop("db must be a 'massFlowTemplate' class object")
            params <- object@params

            while (any(object@samples$aligned == FALSE)) {

              message(paste(length(which(object@samples$aligned == F)), "out of" , nrow(object@samples), "samples left to align."))
              doi_fname <- checkNEXT(object)
              message("Aligning to sample: ", basename(doi_fname),  "... ")
              out <- addDOI(tmp = object@tmp, doi_fname = doi_fname, mz_err = params$mz_err, rt_err = params$rt_err, bins = params$bins)
              object@tmp <- out$tmp
              object@samples[object@samples$filepaths == doi_fname,"aligned"] <- TRUE
              object@data[[doi_fname]] <- out$doi

            }
            message("Peaks were aligned across all samples.")
            return(object)
          }

)