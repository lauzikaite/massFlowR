setMethod("show", signature = "massFlowTemplate", function(object) {
  cat("A \"massFlowTemplate\" object with", nrow(object@samples), " samples")
})

#' Get the the absolute path to the experimental details file, which was used to build the template
#'
#' @param A \code{massFlowTemplate} class object
#'
#' @export
#'
setMethod("filepath", signature = "massFlowTemplate", function(object) object@filepath)

#' Align LC-MS experiment samples using groups of structurally-related peaks
#'
#' @description Merge peak-groups together across samples in the experiment using spectral similarity comparison.
#'
#' @details A sample is aligned to the template, i.e. list of all previously detected peaks, by merging its peaks to the template peaks.
#' Matching peaks are found through \emph(m/z) and \emph(rt) windows. To identify true matches, the overall spectral similarity between a peak-group of the sample and all matching template peak-groups is compared.
#' Spectral similarity is measured by obtaining the cosine of the angle between two 2D vectors, representing each peak-group's \emph(m/z) and \emph(intensity) values. A peak-group of the sample is merged/grouped with the template's peak-group with which the highest cosine was obtained.
#'
#' @param A \code{massFlowTemplate} class object, created by \emph{buildTMP} constructor function.
#'
#' @return A \code{massFlowTemplate} class object with updated \code{tmp}, \code{data} slots.
#' Method also writes intermediate alignment results for every sample in a csv file.
#'
#' @export
#'
setMethod("alignSAMPLES",
          signature = "massFlowTemplate",
          function(object) {

            if (class(object) != "massFlowTemplate") stop("db must be a 'massFlowTemplate' class object")
            params <- object@paramas

            while (any(object@samples$aligned == FALSE)) {

              message(paste(length(which(object@samples$aligned == F)), "out of" , nrow(object@samples), "samples left to align."))
              doi_fname <- checkSAMPLES(object)
              message("Aligning sample: ", basename(doi_fname),  "... ")
              out <- addDOI(tmp = object@tmp, doi_fname = doi_fname, mz_err = params$mz_err, rt_err = params$rt_err, bins = params$bins)
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
                message(paste("Waiting for file to be written for:", doi_fname, "..."))
                Sys.sleep(time = 60)
              }
              doi_present <- file.exists(doi_fname)
            }

            return(doi_fname)
          }
)

