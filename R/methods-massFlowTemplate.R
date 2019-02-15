# show ------------------------------------------------------------------------------------------------------
#' @include classes.R
#' 
#' @rdname massFlowTemplate-class
#' 
#' @param object \code{massFlowTemplate} class object
#' 
#' @export
#' 
setMethod("show", signature = "massFlowTemplate", function(object) {
  cat("A \"massFlowTemplate\" object with",
      nrow(object@samples),
      " samples")
  })

# filepath ------------------------------------------------------------------------------------------------------
#' @include classes.R
#' 
#' @rdname massFlowTemplate-class
#'  
#' @title Obtain absolute path to metadata file of the experiment.
#' 
#' @description Function returs the absolute path to the \code{csv} file used when building the \code{massFlowTemplate} object.
#' File contains study sample names and their acquisition (run) order.
#' 
#' @export
#' 
setMethod("filepath", signature = "massFlowTemplate", function(object) {
  object@filepath
  })

# checkNEXT ------------------------------------------------------------------------------------------------------
#' @title Extract and check the filename of the sample to be processed next
#' 
#' @rdname massFlowTemplate-class
#'
#' @return Method returns filename to be processed next.
#'
#' 
setMethod("checkNEXT",
          signature = "massFlowTemplate",
          function(object) {
            ## extract next-in-line sample and check if it is present (i.e. has been written already)
            doi_fname <-
              object@samples$filepaths[which(object@samples$aligned == FALSE)[1]]
            doi_present <-
              file.exists(doi_fname) # if FALSE, then loop until TRUE
            
            while (!doi_present) {
              if (!file.exists(doi_fname)) {
                message(paste("Waiting for file to be written for:", doi_fname, "..."))
                Sys.sleep(time = 60)
              }
              doi_present <- file.exists(doi_fname)
            }
            return(doi_fname)
          })

# alignPEAKS ------------------------------------------------------------------------------------------------------
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
#' @param object \code{massFlowTemplate} class object, created by \code{buildTMP} constructor function.
#'
#' @return Method updates \code{\link{massFlowTemplate}} class object.
#' Slot @@tmp is updated after each round of sample alignment.
#' Slot @@data stores alignment results for each sample in the experiment.
#' Method writes alignment results to a csv file for each sample in the experiment.
#' 
#' @seealso \code{\link{buildTMP}}, \code{\link{validPEAKS}}, \code{\link{massFlowTemplate-class}}
#'
#' @export
#'
setMethod("alignPEAKS",
          signature = "massFlowTemplate",
          function(object, write_int = TRUE) {
            if (class(object) != "massFlowTemplate") {
              stop("object must be a 'massFlowTemplate' class object")
            }
            params <- object@params
            tmp_fname <-
              gsub(".csv", "_template.csv", object@filepath)
            
            while (any(object@samples$aligned == FALSE)) {
              message(paste(
                length(which(object@samples$aligned == F)),
                "out of" ,
                nrow(object@samples),
                "samples left to align."
              ))
              doi_fname <- checkNEXT(object)
              message("Aligning to sample: ", basename(doi_fname),  "... ")
              out <- addDOI(
                tmp = object@tmp,
                tmp_fname = tmp_fname,
                doi_fname = doi_fname,
                mz_err = params$mz_err,
                rt_err = params$rt_err,
                bins = params$bins,
                write_int = write_int
              )
              object@tmp <- out$tmp
              object@samples[object@samples$filepaths == doi_fname, "aligned"] <-
                TRUE
              object@samples[object@samples$filepaths == doi_fname, "aligned_filepaths"] <-
                gsub(".csv", "_aligned.csv", doi_fname)
              object@data[[doi_fname]] <- out$doi

            }
            ## write updated meta file with filepaths to aligned samples
            write.csv(object@samples, file = gsub(".csv", "_aligned.csv", object@filepath), row.names = F)
            
            message("Peaks were aligned across all samples.")
            return(object)
          })

# validPEAKS ------------------------------------------------------------------------------------------------------
#' @aliases validPEAKS
#'
#' @title Validate aligned peaks and their corresponding peak-groups
#'
#' @param object \code{massFlowTemplate} class object.
#' @param ncores \code{numeric} defining number of cores to use for parallelisation. Default set to 1 for serial implementation.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param min_samples \code{numeric} specifying the minimum percentage of samples in which peak has to be detected in order to be considered (default set to 10 percent). 
#' 
#' @return Method returns validated peak-groups.
#' 
#' @seealso \code{\link{buildTMP}}, \code{\link{alignPEAKS}}, \code{\link{massFlowTemplate-class}}
#'
#' @export
#'
setMethod("validPEAKS",
          signature = "massFlowTemplate",
          function(object, ncores = 1, out_dir, min_samples = 10) {
            
            if (class(object) != "massFlowTemplate") {
              stop("object must be a 'massFlowTemplate' class object")
            }
            if (!is.numeric(ncores)) {
              stop("'ncores' has to be numeric value!")
            }
            if (ncores < 1) {
              stop("'ncores' must be set to 1 (serial performance), or higher!")
            }
            
            message("'ncores' set to ", ncores)
            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)
            if (ncores > 1) {
              bpparam <- BiocParallel::DoparParam()
            } else {
              bpparam <- BiocParallel::SerialParam()
            }
            
            ## extract intensities for every peak-group from every sample
            peakgrs <- unique(object@tmp$peakgr)
            peakgrs <- peakgrs[order(peakgrs)]
            peakgrs_ints <- BiocParallel::bplapply(peakgrs,
                                                   # FUN = extractINT,
                                                   FUN = extractPEAKGR,
                                                   object = object,
                                                   BPPARAM = bpparam)
            saveRDS(peakgrs_ints, file = paste0(out_dir, "/peakgrs_ints.RDS"))
            saveRDS(object, file = paste0(out_dir, "/object.RDS"))
            
            ## validate peak-groups by sample intensity correlation network 
            ## in case validation is run before full alignment of the study, adjust max samples
            samples <- object@samples[which(object@samples$aligned == TRUE),]
            samples_n <- nrow(samples)
            min_samples_n <- ceiling((samples_n * min_samples) / 100)
            peakgrs_split <- BiocParallel::bplapply(
              peakgrs,
              FUN = validPEAKS_paral,
              peakgrs_ints = peakgrs_ints,
              out_dir = out_dir,
              BPPARAM = bpparam,
              min_samples_n = min_samples_n
            )
            saveRDS(peakgrs_split, file = paste0(out_dir, "/peakgrs_split.RDS"))
            
            ## retain only peak-group splits with >=2 peaks
            peakgrs_split_sel <- lapply(peakgrs_split, function(pkg) {
              lapply(pkg, function(pkgs) {
                pkgs_peaks <- length(unique(pkgs$peakid))
                if (pkgs_peaks >= 2) {
                  return(pkgs)
                }
              })
            })
            peakgrs_split_sel <- Filter(length, unlist(peakgrs_split_sel, recursive = F))
            
            ## obtain peakids from the splits
            final_tmp <- lapply(1:length(peakgrs_split_sel), function(n) {
              spkg <- peakgrs_split_sel[[n]]
              data.frame(peakid = unique(spkg$peakid),
                         pcs = n,
                         stringsAsFactors = F)
            })
            final_tmp <- do.call("rbindCLEAN", final_tmp)
            
            ## export intensity measures for each peak and each sample
            final_data <- BiocParallel::bplapply(samples$filename,
                                                 FUN = extractPCS,
                                                 object = object,
                                                 final_tmp = final_tmp
                                                 )
            final_data <- cbind(final_tmp[,c("peakid", "pcs")], 
                                do.call("cbind", final_data))
            
            ## now add remaining centwave output info required for nPYc pipeline
            # final_data <- BiocParallel::bplapply(final_tmp$peakid,
            #                                      FUN = extractPeaks,
            #                                      object = object)
            
            
            
            
            
            
            
            
            parallel::stopCluster(cl)
            
            write.csv(x = final_data,
                      file = gsub(".csv", "_intensity_data.csv", object@filepath),
                      row.names = F)
            write.csv(x = samples,
                      file = gsub(".csv", "_sample_data.csv", object@filepath),
                      row.names = F)
            write.csv(x = final_tmp,
                      file = gsub(".csv", "_peaks_data.csv", object@filepath),
                      row.names = F)
            
            
                                
            message("All peak-groups were succesfully validated.")
            return(object)
          })
