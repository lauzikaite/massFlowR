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
  cat(
    "A \"massFlowTemplate\" object with",
    nrow(object@samples),
    " samples"
  )
})

# setValidity -----------------------------------------------------------------------------------------------------
setValidity("massFlowTemplate", function(object)
  validmassFlowTemplate(object))

# filepath ------------------------------------------------------------------------------------------------------
#' @aliases filepath
#'
#' @title Obtain the absolute path to metadata file of the experiment
#'
#' @description Obtain the absolute path to metadata file of the experiment.
#'
#' @param object \code{massFlowTemplate} class object
#'
#' @seealso \code{\link{massFlowTemplate-class}}
#'
#' @export
#'
setMethod("filepath", signature = "massFlowTemplate", function(object) {
  object@filepath
})

# peaksVALIDATED ------------------------------------------------------------------------------------------------------
#' @aliases peaksVALIDATED
#'
#' @title Check if massFlowTemplate object was validated
#'
#' @description Method returns TRUE if \code{massFlowTemplate} object has been validated via \code{\link{validPEAKS}} method.
#'
#' @param object \code{massFlowTemplate} class object
#'
#' @seealso \code{\link{massFlowTemplate-class}}
#'
#' @export
#'
setMethod("peaksVALIDATED", signature = "massFlowTemplate", function(object) {
  if (nrow(object@valid) > 0) {
    peaks_validated <- TRUE
    if (length(object@peaks) != nrow(object@valid)) {
      peaks_validated <- "Slot 'peaks' doesn't contain a list for every validated peak"
    }
    if (length(object@values) != length(which(object@samples$aligned))) {
      peaks_validated <- "Slot 'values' doesn't contain a list of peak values for every sample"
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
  }
)

# alignPEAKS ------------------------------------------------------------------------------------------------------
#' @aliases alignPEAKS
#'
#' @title Align peaks detected in LC-MS samples using spectral similarity comparison
#'
#' @description Method aligns peaks across samples in LC-MS experiment using spectral similarity comparison.
#' To enable alignment, peaks originating from the same chemical compound were grouped into peak-groups, via function \code{\link{groupPEAKS}}.
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
#' @param ncores \code{numeric} for number of parallel workers to be used. Set 1 for serial implementation. Default set to 2.
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
             ncores = 2,
             write_int = TRUE) {
    if (!validObject(object)) {
      stop(validObject(object))
    }
    if (is.null(out_dir)) {
      stop("'out_dir' is required!")
    }
    if (!dir.exists(out_dir)) {
      stop("incorrect filepath for 'out_dir' provided")
    }

    ## register paral backend
    if (ncores > 1) {
      doParallel::registerDoParallel(cores = ncores)
    } else {
      foreach::registerDoSEQ()
    }
    params <- object@params

    ## filenames for template and final metadata tables
    tmp_fname <- file.path(out_dir, "template.csv")
    aligned_fname <- file.path(out_dir, "aligned.csv")

    while (any(object@samples$aligned == FALSE)) {
      ## check if next peak-group table has been written already
      doi_fname <- checkNEXT(object)
      ## output filenames
      doi_name <- object@samples$filename[object@samples$proc_filepath == doi_fname]
      doi_fname_out <- paste0(file.path(out_dir, doi_name), "_aligned.csv")
      ## check if next peak-group table contains correct columns, if so, read csv
      doi <- checkFILE(file = doi_fname)
      ## status messages
      message(paste(
        length(which(object@samples$aligned == F)),
        "out of",
        nrow(object@samples),
        "samples left to align."
      ))
      message("Aligning to sample: ", doi_name, " ... ")
      ## get most up-to-date template
      tmp <- object@tmp

      doi_to_tmp <- do_alignPEAKS(
        ds = doi,
        tmp = tmp,
        ds_var_name = "peakgr",
        tmp_var_name = "peakgr",
        mz_err = params$mz_err,
        rt_err = params$rt_err,
        bins = params$bins,
        ncores = ncores
      )

      #### ---- update template with aligned doi peak-groups
      ## list to store tmp peakids for every peak in doi
      doi_to_tmp_peakids <- rep(list(NA), nrow(doi))
      doi_to_tmp_peakgr <- rep(list(NA), nrow(doi))
      doi_to_tmp_cos <- rep(list(NA), nrow(doi))

      ## for EVERY peak-group in doi
      for (var in 1:length(doi_to_tmp)) {
        doi_var <- doi_to_tmp[[var]]$ds
        doi_var_peaks <- doi[doi$peakgr == doi_var, ]
        tmp_var <- doi_to_tmp[[var]]$tmp

        if (is.null(tmp_var)) {
          ## use all peaks in the peak-group
          doi_peakids <- doi_var_peaks$peakid
          ## update tmp: add new peaks under new peakgr, set cos to NA
          tmp_peakids <- (max(tmp$peakid) + 1):(max(tmp$peakid) + nrow(doi_var_peaks))
          tmp_var <- max(tmp$peakgr) + 1
          cos <- NA
          utmp <- doi_var_peaks[, c("mz", "rt", "into")]
          utmp$peakid <- tmp_peakids
          utmp$peakgr <- tmp_var
          tmp <- rbind(tmp, utmp)
        } else {
          doi_peakids <- doi_to_tmp[[var]]$mat$target_peakid
          tmp_peakids <- doi_to_tmp[[var]]$mat$peakid
          cos <- doi_to_tmp[[var]]$cos
          ## get mean mz/rt across doi & tmp. Get into values from doi
          new_values <- lapply(1:length(tmp_peakids), function(x) {
            tmp_p <- tmp_peakids[x]
            doi_p <- doi_peakids[x]
            tmp_peak <- tmp[match(tmp_p, tmp$peakid), c("mz", "rt", "into")]
            doi_peak <- doi[match(doi_p, doi$peakid), c("mz", "rt", "into")]
            rt_m <- mean(c(doi_peak$rt, tmp_peak$rt))
            mz_m <- mean(c(doi_peak$mz, tmp_peak$mz))
            into <- doi_peak$into[which.max(doi_peak$into)]
            data.frame(mz_m, rt_m, into)
          })
          tmp[match(tmp_peakids, tmp$peakid), c("mz", "rt", "into")] <- do.call("rbind", new_values)
          ## if any unmatched peaks are left, add them as new
          doi_var_peaks <- doi[which(doi$peakgr == doi_var), ]
          doi_var_peaks_un <- doi_var_peaks[!doi_var_peaks$peakid %in% doi_peakids, ]
          if (nrow(doi_var_peaks_un) > 0) {
            tmp_peakids_new <- (max(tmp$peakid) + 1):(max(tmp$peakid) + nrow(doi_var_peaks_un))
            utmp <- doi_var_peaks_un[, c("mz", "rt", "into")]
            utmp$peakid <- tmp_peakids_new
            utmp$peakgr <- tmp_var
            tmp <- rbind(tmp, utmp)
            tmp_peakids <- c(tmp_peakids, tmp_peakids_new)
            doi_peakids <- c(doi_peakids, doi_var_peaks_un$peakid)
          }
        }
        ## update doi-to-tmp matching info
        doi_to_tmp_peakids[unlist(doi_peakids)] <- tmp_peakids
        doi_to_tmp_peakgr[unlist(doi_peakids)] <- tmp_var
        doi_to_tmp_cos[unlist(doi_peakids)] <- cos
      }
      doi$tmp_peakid <- unlist(doi_to_tmp_peakids)
      doi$tmp_peakgr <- unlist(doi_to_tmp_peakgr)
      doi$cos <- unlist(doi_to_tmp_cos)

      object@tmp <- tmp
      object@samples[object@samples$proc_filepath == doi_fname, "aligned"] <-
        TRUE
      object@samples[object@samples$proc_filepath == doi_fname, "aligned_filepath"] <-
        doi_fname_out
      object@data[[doi_name]] <- doi

      ## (1) write aligned doi table
      write.csv(
        doi,
        file = doi_fname_out,
        quote = T,
        row.names = F
      )
      ## (2) write intermediate template generated after this doi
      if (write_int == T) {
        write.csv(
          tmp,
          file = gsub("aligned.csv", "tmp.csv", doi_fname_out),
          quote = T,
          row.names = F
        )
      }
      ## (3) overwrite template file
      write.csv(tmp,
        file = tmp_fname,
        quote = T,
        row.names = F
      )
    }
    object@history[length(object@history) + 1] <- "alignPEAKS"
    ## write updated meta file with filepaths to aligned samples
    write.csv(object@samples,
      aligned_fname,
      quote = T,
      row.names = F
    )

    if (ncores > 1) {
      foreach::registerDoSEQ()
    }
    message("Peaks were aligned across all samples.")
    return(object)
  }
)

# validPEAKS ------------------------------------------------------------------------------------------------------
#' @aliases validPEAKS
#'
#' @title Validate aligned peaks and their corresponding peak-groups
#'
#' @param object \code{massFlowTemplate} class object.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param min_samples_n \code{numeric} specifying the minimum percentage of samples in which peak has to be detected in order to be considered (default set to 10 percent).
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
             min_samples_n = 3,
             cor_thr = 0.75,
             ncores = 2) {
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
      doParallel::registerDoParallel(cores = ncores)
    } else {
      foreach::registerDoSEQ()
    }

    ## extract intensities for every peak-group from every sample
    peakgrs <- unique(object@tmp$peakgr)
    peakgrs <- peakgrs[order(peakgrs)]
    peakgrs_ints <- foreach::foreach(
      pkg = peakgrs,
      .inorder = TRUE,
      .export = c("extractPEAKGR")
    ) %dopar% (extractPEAKGR(
      pkg = pkg,
      object = object
    ))
    saveRDS(peakgrs_ints, file = paste0(out_dir, "/peakgrs_ints.RDS"))
    saveRDS(object, file = paste0(out_dir, "/object.RDS"))

    ## validate peak-groups by sample intensity correlation network
    ## in case validation is run before full alignment of the study, adjust sample n
    samples <-
      object@samples[which(object@samples$aligned == TRUE), ]
    samples_n <- nrow(samples)
    if (min_samples_n < 3) {
      min_samples_n <- 3
    }
    if (min_samples_n > samples_n) {
      stop(
        "object has ", samples_n, " samples",
        "\n minimum 3 samples are required for validation."
      )
    }
    peakgrs_split <- foreach::foreach(
      pkg = seq(length(peakgrs)),
      .inorder = TRUE,
      .export = c("validPEAKGR")
    ) %dopar% (
      validPEAKGR(
        pkg = pkg,
        pkg_ints = peakgrs_ints[[pkg]],
        out_dir = out_dir,
        min_samples_n = min_samples_n,
        cor_thr = cor_thr
      )
    )
    saveRDS(peakgrs_split, file = paste0(out_dir, "/peakgrs_split.RDS"))

    ## retain only communities with > 1 peak
    final_tmp <-
      extractCOMMUNITIES(peakgrs_split = peakgrs_split)

    ## export centWave measures for final template's peaks in each sample, listed by sample
    peaks_vals_samples <-
      foreach::foreach(
        sn = 1:nrow(samples),
        .inorder = TRUE,
        .export = c("exportSAMPLE")
      ) %dopar% (exportSAMPLE(
        sdata = object@data[[sn]],
        final_tmp = final_tmp
      ))

    ## get centWave measures for each peak and each sample, listed by peakid
    peaks_vals_peakids <-
      foreach::foreach(
        peakid = final_tmp$peakid,
        .inorder = TRUE,
        .export = c("exportPEAK")
      ) %dopar% (exportPEAK(
        peakid = peakid,
        peaks_vals_samples = peaks_vals_samples
      ))
    names(peaks_vals_peakids) <- final_tmp$peakid

    ## get median centWave measures for each peak and each sample, listed by peakid
    peaks_medians_peakids <-
      foreach::foreach(
        n = 1:length(peaks_vals_peakids),
        .inorder = TRUE,
        .export = c("getPEAKmedians")
      ) %dopar% (getPEAKmedians(peak_n = peaks_vals_peakids[[n]]))
    peaks_medians_peakids <-
      do.call("rbindCLEAN", peaks_medians_peakids)
    peaks_medians_peakids$peakid <- final_tmp$peakid
    peaks_medians_peakids$pcs <-
      final_tmp$pcs[match(peaks_medians_peakids$peakid, final_tmp$peakid)]

    ## merge intensity and centWave values into one dataframe
    ## replace NA values with 0, extract for each sample
    peaks_vals <-
      lapply(1:length(peaks_vals_samples), function(sn) {
        sname <- samples$filename[sn]
        sdata <- peaks_vals_samples[[sn]]
        sdata_out <-
          as.data.frame(ifelse(is.na(sdata$into), 0, sdata$into),
            row.names = NULL
          )
        sdata_out <- setNames(sdata_out, sname)
        return(sdata_out)
      })
    peaks_vals <- do.call("cbind", peaks_vals)
    peaks_vals <- cbind(peaks_medians_peakids, peaks_vals)

    object@values <- peaks_vals_samples
    object@peaks <- peaks_vals_peakids
    object@valid <- peaks_medians_peakids
    object@history[length(object@history) + 1] <- "validPEAKS"

    #### ---- data output
    write.csv(
      x = peaks_vals,
      file = file.path(out_dir, "intensity_data.csv"),
      row.names = F
    )
    write.csv(
      x = samples,
      file = file.path(out_dir, "sample_data.csv"),
      row.names = F
    )
    write.csv(
      x = final_tmp,
      file = file.path(out_dir, "peaks_data.csv"),
      row.names = F
    )
    if (ncores > 1) {
      foreach::registerDoSEQ()
    }
    message("All peak-groups were succesfully validated.")
    return(object)
  }
)

# fillPEAKS ------------------------------------------------------------------------------------------------------
#' @aliases fillPEAKS
#'
#' @title Fill peaks
#'
#' @param object \code{massFlowTemplate} class object.
#' @param fill_value \code{character} specifying which intensity value should be filled and returned, default set to 'into'.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param ncores \code{numeric} defining number of cores to use for parallelisation. Default set to 2 for parallel implementation.
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
             fill_value = "into",
             out_dir = NULL,
             ncores = 2) {
    if (!validObject(object)) {
      stop(validObject(object))
    }
    if (!peaksVALIDATED(object) | !"validPEAKS" %in% object@history) {
      stop(
        "'object' must be a validated 'massFlowTemplate' class object. \n Run validPEAKS() first."
      )
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
    ## check if raw files are valid before proceeding
    snames_valid <- unlist(lapply(object@samples$raw_filepath, FUN = validFILE))
    if (length(which(snames_valid) == TRUE) != length(object@samples$raw_filepath)) {
      stop(
        "raw filepaths in massFlowTemplate object are incorrect: ",
        snames_valid[which(snames_valid != TRUE)]
      )
    }
    ## create paralel backend initiation function to be used in cluster re-initiation
    if (ncores > 1) {
      registerBACKEND <- function(ncores) {
        doParallel::registerDoParallel(cores = ncores)
      }
    } else {
      ## if serial implementation is required, register sequential backend (to avoid warning message)
      registerBACKEND <- function(ncores) {
        foreach::registerDoSEQ()
      }
    }
    ## register backend for modelPEAKS and extractMISS functions
    registerBACKEND(ncores = ncores)

    #### ---- model mz/rt windows for missing peaks
    ## assuming that every peak will have at least one sample to fill,
    ## iterate over every peakid in the validated peak-table
    peakids <- object@valid$peakid
    peaks_mod <- foreach::foreach(
      p = peakids,
      .inorder = TRUE,
      .export = c("modelPEAKS")
    ) %dopar% (modelPEAKS(
      p = p,
      vars = c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax"),
      object = object
    ))

    #### ---- fill missing peaks
    ## build list of mz/rt regions to fill, sample-wise
    ## extract only the peaks that have to be filled
    peaks_miss <-
      foreach::foreach(
        s = 1:length(object@values),
        .inorder = TRUE,
        .export = c("extractMISS")
      ) %dopar% (
        extractMISS(
          s = s,
          values = object@values,
          peaks_mod = peaks_mod,
          valid = object@valid
        )
      )
    ## updates centWave values for every sample with newly integrated into and maxo values
    samples_to_fill <- 1:length(object@values)

    ## create named vector to store processing result for every sample
    ## this list will be reused for failed samples as well
    nsamples <- length(samples_to_fill)
    result <- vector("list", nsamples)

    while (length(samples_to_fill) > 0) {

      ## update nsamples
      nsamples <- length(samples_to_fill)

      ## get number of processes across which samples will be divided
      nproc <- ceiling(nsamples / ncores)

      ## in each new process, start a new cluster with selected number of cores
      for (iproc in seq(0, (nproc - 1), by = 1)) {

        ## re-intialise cluster for this set of files (helps reduce memory usage)
        registerBACKEND(ncores = ncores)

        ## extract files to be processed in this process
        samples_proc_first <- 1 + (iproc * ncores)
        samples_proc_last <- ncores + (iproc * ncores)
        samples_proc_last <- ifelse(samples_proc_last <= nsamples, samples_proc_last, nsamples)
        samples_proc <- samples_to_fill[samples_proc_first:samples_proc_last]
        message(length(samples_to_fill), " samples left to fill.")
        message(
          "Filling next ", length(samples_proc), " samples: \n",
          paste0(object@samples$filename[samples_proc], collapse = "\n"),
          "\n ..."
        )
        result[samples_proc] <- foreach::foreach(
          s = samples_proc,
          .inorder = TRUE,
          .errorhandling = "pass",
          .export = c("fillSAMPLE")
        ) %dopar% (tryCatch({
          fillSAMPLE(
            s = s,
            sname = object@samples$raw_filepath[s],
            sdata = object@values[[s]],
            values = peaks_miss[[s]]
          )
        }, error = function(err) {
          return(list(
            nsample = s,
            status = "FAILED",
            error = err$message
          ))
        }))
      }

      ## process results
      # identify files that failed and update list of files
      result_status <- sapply(result, "[[", "status")
      samples_to_fill <- which(result_status == "FAILED")

      if (length(samples_to_fill) > 0) {
        message(" \n", length(samples_to_fill), " samples failed. Filling failed samples ...")
        message("Generated error messages for failed samples:\n")
        message(paste0(sapply(result[samples_to_fill], "[[", "error")))
      }
    }

    ## update peak medians with median into/maxo values after filling
    peaks_vals <-
      lapply(1:length(object@values), function(sn) {
        sname <- object@samples$filename[sn]
        sdata <- result[[sn]]
        sdata_out <-
          as.data.frame(sdata[, match(fill_value, colnames(sdata))],
            row.names = NULL
          )
        sdata_out <- setNames(sdata_out, sname)
        return(sdata_out)
      })

    peaks_vals <- do.call("cbind", peaks_vals)
    peaks_median_vals <- apply(peaks_vals, 1, function(x) {
      median(x[which(x > 0)])
    })
    tmp <- object@valid
    tmp[, "into"] <- peaks_median_vals
    peaks_vals <- cbind(tmp, peaks_vals)

    ## replace object's values
    object@values <- result
    object@valid <- tmp
    object@history[length(object@history) + 1] <- "fillPEAKS"

    #### ---- data output
    write.csv(
      x = peaks_vals,
      file = file.path(out_dir, "filled_intensity_data.csv"),
      row.names = F
    )
    write.csv(
      x = tmp,
      file = file.path(out_dir, "final_peaks_data.csv"),
      row.names = F
    )
    if (validObject(object)) {
      message("All peak-groups were succesfully filled")
      return(object)
    }
  }
)

# adjustBATCH -------------------------------------------------------------
#' @aliases adjustBATCH
#'
#' @title adjustBATCH
#'
#' @description Development mode.
#' Method adjusts retention time between two analytical batches using a set of features for which regions of integration (ROI) in \emph{mz}, \emph{rt} and \emph{into} domains are specified by the user.
#'
#' @param object \code{massFlowTemplate} class object.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param batch_end \code{character} with filename of the last sample in batch No 1.
#' @param batch_start \code{character} with filename of the first sample in batch No 2.
#' @param batch_next_metadata \code{character} with path to metadata file for batch No 2.
#' @param batch_end_roi \code{character} with path to csv file with regions of integrations for the last sample in batch No 1.
#' @param batch_start_roi \code{character} with path to csv file with regions of integrations for the first sample in batch No 2.
#'
#' @return Method returns \code{massFlowTemplate} class object with updated template, such that \emph{rt} of the features are adjusted according to the observed batch-change effect.
#'
#' @export
#'
setMethod("adjustBATCH",
  signature = "massFlowTemplate",
  function(object,
             out_dir = NULL,
             batch_end = NULL,
             batch_start = NULL,
             batch_next_metadata = NULL,
             batch_end_roi = NULL,
             batch_start_roi = NULL) {
    if (!validObject(object)) {
      stop(validObject(object))
    }
    # if (!peaksVALIDATED(object)) {
    #   stop("object must be validated first")
    # }
    if (any(is.null(out_dir) |
      is.null(batch_end) |
      is.null(batch_start) |
      is.null(batch_next_metadata) |
      is.null(batch_end_roi) |
      is.null(batch_start_roi))) {
      stop("check your arguments")
    }
    if (!dir.exists(out_dir)) {
      stop("incorrect filepath for 'out_dir' provided")
    }
    params <- object@params
    ## load and check provided reference compounds tables
    roi_end <-
      read.csv(batch_end_roi,
        header = TRUE,
        stringsAsFactors = FALSE
      )
    roi_start <-
      read.csv(batch_start_roi,
        header = TRUE,
        stringsAsFactors = FALSE
      )
    if (nrow(roi_end) != nrow(roi_start)) {
      stop("provided reference compounds integration tables are inconsistent")
    }
    ## load metadata of the next batch
    next_samples <- read.csv(batch_next_metadata, header = TRUE, stringsAsFactors = FALSE)

    ## final template with all valid peaks
    # val_pks <- object@valid # validated peaks
    val_pks <- object@tmp # fix for simulated data

    ## load peak tables for the end of the current batch and start of the next batch samples
    data_end <-
      object@data[[which(object@samples$filename == batch_end)]] # batch end sample - user selected file
    # data_end <- object@data[[94]] # lazy fix for devset
    data_start <-
      read.csv(next_samples$proc_filepath[which(next_samples$filename == batch_start)],
        header = TRUE,
        stringsAsFactors = FALSE
      ) # starting sample of the next batch
    rt_difs <- lapply(
      1:nrow(roi_end),
      FUN = checkCOMPOUND,
      val_pks = val_pks,
      roi_end = roi_end,
      roi_start = roi_start,
      data_end = data_end,
      data_start = data_start
    )
    if (length(unlist(rt_difs)) <= 3) {
      ans <- 0
      while (ans < 1) {
        ans <- readline(
          paste0(
            "Peaks corresponding to provided reference compounds: ",
            length(unlist(rt_difs)),
            "\n",
            "Do you wish to proceed with batch adjustion? Enter Y/N "
          )
        )
        ## catch if input is N/n
        ans <- ifelse((grepl("N", ans) | grepl("n", ans)),
          2, 1
        )
        if (ans == 2) {
          stop("method was stopped.")
        }
      }
    }
    ## prepare template for next batch:
    ## (1) extract validated peaks
    ## (2) adjust their rt values
    tmp <- object@tmp
    tmp <- tmp[match(val_pks$peakid, tmp$peakid), ]
    ## version A - median of rt differences
    rt_dif <- median(unlist(rt_difs))

    ## version B - model rt change
    graphics::plot(roi_end$rt, unlist(lapply(rt_difs, function(x) {
      ifelse(is.null(x), NA, x)
    })))

    tmp$rt <- tmp$rt + rt_dif

    ## create new object for the next batch, which will replace current object
    next_samples[, "aligned"] <- FALSE
    next_samples[, "aligned_filepath"] <- NA
    next_object <- new(
      "massFlowTemplate",
      filepath = batch_next_metadata,
      samples = next_samples,
      tmp = tmp,
      params = list(
        mz_err = params$mz_err,
        rt_err = params$rt_err,
        bins = params$bins
      )
    )
    if (validmassFlowTemplate(next_object) != TRUE) {
      stop(validmassFlowTemplate(next_object))
    }
    return(next_object)
  }
)
