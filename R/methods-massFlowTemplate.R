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


# setValidity -----------------------------------------------------------------------------------------------------
setValidity("massFlowTemplate", function(object) validmassFlowTemplate(object))


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


# peaksVALIDATED ------------------------------------------------------------------------------------------------------
#' @include classes.R
#' 
#' @rdname massFlowTemplate-class
#'
#' @title Check if massFlowTemplate object was validated using validPEAKS.
#'
#' @export
#'
setMethod("peaksVALIDATED", signature = "massFlowTemplate", function(object) {
  if (nrow(object@valid) > 0) {
    peaks_validated <- TRUE
    if (length(object@peaks) != nrow(object@valid)) {
      msg <- c(msg, paste0("Slot 'peaks' doesn't contain a list for every validated peak"))
    }
    if (length(object@values) != length(which(object@samples$aligned))) {
      msg <- c(msg, paste0("Slot 'values' doesn't contain a list of peak values for every sample"))
    }
  } else {
    peaks_validated <- FALSE
  }
  return(peaks_validated)
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
              object@samples$proc_filepath[which(object@samples$aligned == FALSE)[1]]
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
#' @param out_dir \code{character} specifying desired directory for output.
#' @param write_int \code{logical} specifying whether a peak table with alignment results should be saved for every sample.
#' If TRUE, csv files will be written in the out_dir directory.
#' Default set to TRUE.
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
          function(object,
                   out_dir = NULL,
                   write_int = TRUE
                   ) {
            if (!validObject(object)) {
              stop(validObject(object))
            }
            if (is.null(out_dir)) {
              stop("'out_dir' is required!")
            }
            if (!dir.exists(out_dir)) {
              stop("incorrect filepath for 'out_dir' provided")
            }
           
            params <- object@params
            
            ## generated template and all intermediate alignment peak tables
            tmp_fname <- file.path(out_dir, "template.csv")
            aligned_fname <- file.path(out_dir, "aligned.csv")
            
            while (any(object@samples$aligned == FALSE)) {
              message(paste(
                length(which(object@samples$aligned == F)),
                "out of" ,
                nrow(object@samples),
                "samples left to align."
              ))
              doi_fname <- checkNEXT(object)
              doi_name <- object@samples$filename[object@samples$proc_filepath == doi_fname]
              doi_fname_out <- paste0(file.path(out_dir, doi_name), "_aligned.csv")
              message("Aligning to sample: ", doi_name,  " ... ")
              out <- addDOI(
                tmp = object@tmp,
                tmp_fname = tmp_fname,
                doi_fname = doi_fname,
                doi_fname_out = doi_fname_out,
                mz_err = params$mz_err,
                rt_err = params$rt_err,
                bins = params$bins,
                write_int = write_int
              )
              object@tmp <- out$tmp
              object@samples[object@samples$proc_filepath == doi_fname, "aligned"] <-
                TRUE
              object@samples[object@samples$proc_filepath == doi_fname, "aligned_filepath"] <-
                doi_fname_out
              object@data[[doi_name]] <- out$doi

            }
            ## write updated meta file with filepaths to aligned samples
            write.csv(
              object@samples,
              aligned_fname,
              quote = T,
              row.names = F
            )
            message("Peaks were aligned across all samples.")
            return(object)
          })

# validPEAKS ------------------------------------------------------------------------------------------------------
#' @aliases validPEAKS
#'
#' @title Validate aligned peaks and their corresponding peak-groups
#'
#' @param object \code{massFlowTemplate} class object.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param min_samples \code{numeric} specifying the minimum percentage of samples in which peak has to be detected in order to be considered (default set to 10 percent). 
#' @param cor_thr \code{numeric} defining Pearson correlation coefficient threshold for inter-sample correlation between peaks (default set to 0.75). 
#' @param ncores \code{numeric} defining number of cores to use for parallelisation. Default set to 1 for serial implementation.
#' 
#' @return Method returns validated peak-groups.
#' 
#' @seealso \code{\link{buildTMP}}, \code{\link{alignPEAKS}}, \code{\link{massFlowTemplate-class}}
#'
#' @export
#'
setMethod("validPEAKS",
          signature = "massFlowTemplate",
          function(object,
                   out_dir = NULL,
                   min_samples = 10,
                   cor_thr = 0.75,
                   ncores = 2
                   ) {
            if (!validObject(object)) {
              stop(validObject(object))
            }
            if (is.null(out_dir)) {
              stop("'out_dir' is required!")
            }
            if (!dir.exists(out_dir)) {
              stop("incorrect filepath for 'out_dir' provided")
            }
            if (ncores < 1 | !is.numeric(ncores)) {
              warning("'ncores' was not correctly set. Switching to ncores = 1 (serial performance)")
            }
            ## register paral backend
            if (ncores > 1) {
              cl <- parallel::makeCluster(ncores)
              doParallel::registerDoParallel(cl)
            } else {
              foreach::registerDoSEQ()
            }
              
            ## extract intensities for every peak-group from every sample
            peakgrs <- unique(object@tmp$peakgr)
            peakgrs <- peakgrs[order(peakgrs)]
            peakgrs_ints <- foreach::foreach(
              pkg = peakgrs,
              .inorder = TRUE
            ) %dopar% (massFlowR:::extractPEAKGR(
              pkg = pkg,
              object = object
            ))
            saveRDS(peakgrs_ints, file = paste0(out_dir, "/peakgrs_ints.RDS"))
            saveRDS(object, file = paste0(out_dir, "/object.RDS"))
            
            ## validate peak-groups by sample intensity correlation network 
            ## in case validation is run before full alignment of the study, adjust sample n
            samples <- object@samples[which(object@samples$aligned == TRUE),]
            samples_n <- nrow(samples)
            min_samples_n <- ceiling((samples_n * min_samples) / 100)
            if (min_samples_n < 3) {
              min_samples_n <- 3
            }
            peakgrs_split <- foreach::foreach(
              pkg = peakgrs,
              .inorder = TRUE
            ) %dopar% (massFlowR:::validPEAKGR(
              pkg = pkg,
              peakgrs_ints = peakgrs_ints,
              out_dir = out_dir,
              min_samples_n = min_samples_n,
              cor_thr = cor_thr
            ))
            saveRDS(peakgrs_split, file = paste0(out_dir, "/peakgrs_split.RDS"))
            
            ## retain only communities with > 1 peak
            final_tmp <- extractCOMMUNITIES(peakgrs_split = peakgrs_split)
            
            ## export centWave measures for final template's peaks in each sample, listed by sample
            peaks_vals_samples <- foreach::foreach(
              sn = 1:nrow(samples),
              .inorder = TRUE
            ) %dopar% (massFlowR:::exportSAMPLE(
              sdata = object@data[[sn]],
              final_tmp = final_tmp
            ))
  
            ## get centWave measures for each peak and each sample, listed by peakid
            peaks_vals_peakids <- foreach::foreach(
              peakid = final_tmp$peakid,
              .inorder = TRUE
            ) %dopar% (massFlowR:::exportPEAK(
                peakid = peakid,
                peaks_vals_samples = peaks_vals_samples
              ))
            names(peaks_vals_peakids) <- final_tmp$peakid

            ## get median centWave measures for each peak and each sample, listed by peakid
            peaks_medians_peakids <- foreach::foreach(
              n = 1:length(peaks_vals_peakids),
              .inorder = TRUE
            ) %dopar% (massFlowR:::getPEAKmedians(
              peak_n = peaks_vals_peakids[[n]]
              ))
            peaks_medians_peakids <- do.call("rbindCLEAN", peaks_medians_peakids)
            peaks_medians_peakids$peakid <- final_tmp$peakid
            peaks_medians_peakids$pcs <-
              final_tmp$pcs[match(peaks_medians_peakids$peakid, final_tmp$peakid)]
          
            ## merge intensity and centWave values into one dataframe
            ## replace NA values with 0, extract for each sample
            peaks_vals <- lapply(1:length(peaks_vals_samples), function(sn) {
              sname <- samples$filename[sn]
              sdata <- peaks_vals_samples[[sn]]
              sdata_out <-
                as.data.frame(ifelse(is.na(sdata$into), 0, sdata$into),
                              row.names = NULL)
              sdata_out <- setNames(sdata_out, sname)
              return(sdata_out)
            })
            peaks_vals <- do.call("cbind", peaks_vals)
            peaks_vals <- cbind(peaks_medians_peakids, peaks_vals)

            parallel::stopCluster(cl)
            
            object@values <- peaks_vals_samples
            object@peaks <- peaks_vals_peakids
            object@valid <- final_tmp
            
            ####---- data output
            write.csv(x = peaks_vals,
                      file = file.path(out_dir, "intensity_data.csv"),
                      row.names = F)
            write.csv(x = samples,
                      file = file.path(out_dir, "sample_data.csv"),
                      row.names = F)
            write.csv(x = final_tmp,
                      file = file.path(out_dir, "peaks_data.csv"),
                      row.names = F)
               
            message("All peak-groups were succesfully validated.")
            return(object)
          })


# fillPEAKS ------------------------------------------------------------------------------------------------------
#' @aliases fillPEAKS
#'
#' @title Fill peaks 
#'
#' @param object \code{massFlowTemplate} class object.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param min_samples \code{numeric} specifying the minimum percentage of samples in which peak has to be detected in order to be considered (default set to 10 percent). 
#' @param cor_thr \code{numeric} defining Pearson correlation coefficient threshold for inter-sample correlation between peaks (default set to 0.75). 
#' @param ncores \code{numeric} defining number of cores to use for parallelisation. Default set to 1 for serial implementation.
#' 
#' @return Method returns peak table with filled intensities for missed peaks.
#' 
#' @seealso \code{\link{loadALIGNED}}, \code{\link{alignPEAKS}}, \code{\link{validPEAKS}}, \code{\link{massFlowTemplate-class}}
#'
#' @export
#'
setMethod("fillPEAKS",
          signature = "massFlowTemplate",
          function(object,
                   out_dir = NULL,
                   ncores = 1
          ) {
            if (!validObject(object)) {
              stop(validObject(object))
            }
            if (!peaksVALIDATED(object)) {
              stop("'Object' must be a validated 'massFlowTemplate' class object. \n Run validPEAKS() first.")
            }
            if (is.null(out_dir)) {
              stop("'out_dir' is required!")
            }
            if (!dir.exists(out_dir)) {
              stop("incorrect filepath for 'out_dir' provided")
            }
            if (ncores < 1 | !is.numeric(ncores)) {
              warning("'ncores' was not correctly set. Switching to ncores = 1 (serial performance)")
            }
            ## register paralel backend
            if (ncores > 1) {
              cl <- parallel::makeCluster(ncores)
              doParallel::registerDoParallel(cl)
            } else {
              foreach::registerDoSEQ()
            }
          
          
            ####---- model mz/rt windows for missing peaks
            ## assuming that every peak will have at least one sample to fill,
            ## iterate over every peakid in the validated peak-table
            peakids <- object@valid$peakid
            peaks_mod <- foreach::foreach(
              p = peakids,
              .inorder = TRUE
            ) %dopar% modelPEAKS(
              p = p,
              vars = c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax"),
              object = object
            )
            
            ####---- fill missing peaks
            ## build list of mz/rt regions to fill, sample-wise
            ## extract only the peaks that have to be filled
            peaks_miss <- foreach::foreach(
              s = 1:length(object@values),
              .inorder = TRUE
            ) %dopar% extractMISS(
              s = s,
              values = object@values,
              peaks_mod = peaks_mod,
              valid = object@valid
            )
            
            peaks_fill <- foreach::foreach(
              s = 1:length(object@values),
              .inorder = TRUE
            ) %dopar% fillPEAKS_sample(
              s = s,
              values = peaks_miss
            )
            
            
            
            
            
            if (validObject(object)) {
              return(object)
            }
          })