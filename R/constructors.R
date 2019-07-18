# buildTMP --------------------------------------------------------------------------------------------------------
#' @title Build a sample alignment and annotation template
#'
#' @description  Functions builds a \code{massFlowTemplate} class object, which stores study sample information.
#'
#' @param file \code{character} for absolute path to the csv file, specifying samples filenames and their acquisition order.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param mz_err \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.05).
#'
#' @return A \code{massFlowTemplate} class object.
#'
#' @export
buildTMP <-
  function(file = NULL,
           out_dir = NULL,
           mz_err = 0.01,
           rt_err = 2,
           bins = 0.05
           ) {
    if (is.null(file)) {
      stop("'file' is required")
    }
    if (!file.exists(file)) {
      stop("incorrect filepath for 'file' provided")
    }
    if (is.null(out_dir)) {
      stop("'out_dir' is required")
    }
    if (!dir.exists(out_dir)) {
      stop("incorrect filepath for 'out_dir' provided")
    }
    samples <- read.csv(file, header = T, stringsAsFactors = F)
    samples[, "aligned"] <- FALSE
    samples[, "aligned_filepath"] <- NA
    object <- new(
      "massFlowTemplate",
      filepath = file,
      samples = samples,
      params = list(
        mz_err = mz_err,
        rt_err = rt_err,
        bins = bins
      )
    )
    if (validmassFlowTemplate(object) != TRUE) {
      stop(validmassFlowTemplate(object))
    }
    doi_first <- min(object@samples$run_order)
    doi_fname <-
      object@samples$proc_filepath[object@samples$run_order == doi_first]
    doi_name <- object@samples$filename[object@samples$run_order == doi_first]
    ## get filename for the file to be written in the selected directory
    doi_fname_out <- paste0(file.path(out_dir, doi_name), "_aligned.csv")
    message(paste("Building template using sample:", doi_name, " ..."))
    ## write 1st sample in the standard output format
    doi <- checkFILE(file = doi_fname)
    doi[, c("tmp_peakid", "tmp_peakgr")] <- doi[, c("peakid", "peakgr")]
    doi[,c("cos")] <- NA
    write.csv(
      doi,
      file = doi_fname_out,
      quote = T,
      row.names = F
    ) 
    ## build template from 1st sample
    tmp <-
      doi[, c("peakid",
              "mz",
              "rt",
              "into",
              "peakgr")]

    object@tmp <- tmp
    object@samples[object@samples$run_order == doi_first, "aligned"] <-
      TRUE
    ## first aligned sample is written in the defined directory
    object@samples[object@samples$run_order == doi_first, "aligned_filepath"] <-
      doi_fname_out
    object@data[[doi_name]] <-
      doi
    return(object)
  }


# loadALIGNED -----------------------------------------------------------------------------------------------------
#' @title Build sample alignment and annotation template using already aligned samples
#'
#' @description Function handles the construction of \code{massFlowTemplate} class object from already aligned samples.
#' Function enables user to continue an interrupted peak alignment process, facilitated by \code{alignPEAKS} function.
#'
#' @details Arguments are identical to the ones used by \code{\link{buildTMP}} constructor function.
#'
#' @param file A \code{character} with path to the csv file, specifying samples filenames and their acquisition order.
#' @param template A \code{character} with path to the csv file with the latest template obtained by \code{alignPEAKS} function.
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 2 (sec).
#' @param bins A \code{numeric} defying step size used in component's spectra binning and vector generation. Step size represents MZ dimension (default set to 0.05).
#'
#' @return A \code{massFlowTemplate} class object.
#'
#' @seealso \code{\link{massFlowTemplate}} class.
#'
#' @export
#'
loadALIGNED <-
  function(file = NULL,
           template = NULL,
           mz_err = 0.01,
           rt_err = 2,
           bins = 0.05) {
  
    if (is.null(file)) {
      stop("Input 'file' is required")
    }
    if (!file.exists(file)) {
      stop("Incorrect filepath for 'file' provided")
    }
    req_cnames <- c("filename",
                    "run_order",
                    "raw_filepath",
                    "proc_filepath",
                    "aligned_filepath")
    samples <- read.csv(file, header = T, stringsAsFactors = F)
    if (any(!req_cnames %in% names(samples))) {
      stop("'files' table must contain columns: ", paste0(req_cnames, collapse = ", "))
    }
    if (any(!file.exists(samples$proc_filepath))) {
      stop("Column 'proc_filepath' contain incorrect file paths: ",
           samples$proc_filepath[!file.exists(samples$proc_filepath)])
    }
    if (any(!file.exists(samples$aligned_filepath))) {
      warning(
        "Column 'aligned_filepath' contain incorrect file paths: ",
        samples$aligned_filepath[!file.exists(samples$aligned_filepath)],
        "\n ",
        "Only correct 'aligned_filepath' will be loaded."
      )
      ans <- 0
      while (ans < 1) {
        ans <- readline("Continue? Enter Y/N ")
        ## catch if input is N/n
        ans <- ifelse((grepl("N", ans) | grepl("n", ans)),
                      2, 1)
        if (ans == 2) {
          stop("loading was stopped.")
        }
      }
    }
    ## load provided template file
    if (!file.exists(template)) {
      stop("template file is not available: ", template)
    }
    ## extract only already aligned samples from the provided file list
    samples_aligned <- which(file.exists(samples$aligned_filepath))
    tmp <- read.csv(template, header = T, stringsAsFactors = F)
    
    object <- new("massFlowTemplate")
    object@filepath <- file
    object@samples <- samples
    object@samples[, "aligned"] <- FALSE
    object@samples[samples_aligned, "aligned"] <- TRUE
    object@params <-
      list(
        mz_err = mz_err,
        rt_err = rt_err,
        bins = bins
      )
    object@tmp <- tmp
    
    ## load aligned samples datasets
    data <-
      lapply(samples$aligned_filepath[samples_aligned], function(doi_fname) {
        doi <- read.csv(doi_fname,
                        header = T,
                        stringsAsFactors = F)
      })
    names(data) <- object@samples[samples_aligned, "proc_filepath"]
    object@data <- data
    message("A 'massFlowTemplate' object was succesfully built with aligned samples.")
    return(object)
  }


# buildANNO ---------------------------------------------------------------
#' @title Build a \code{massFlowAnno} class object.
#'
#' @param ds_file \code{character} for absolute path to the csv file with final peak table obtained with \code{\link{fillPEAKS}}.
#' @param meta_file \code{character} for absolute path to the csv file, specifying samples filenames and their acquisition order.
#' @param out_dir \code{character} specifying desired directory for output.
#' 
#' @return A \code{massFlowAnno} class object.
#' 
#' @export
#'
buildANNO <- function(ds_file = NULL,
                      meta_file = NULL,
                      out_dir = NULL
                      ) {
  
  if (is.null(ds_file) | is.null(meta_file)) {
    stop("'ds_file' and  'meta_file' are required")
  }
  if (!file.exists(ds_file)) {
    stop("incorrect filepath for 'ds_file' provided")
  }
  if (!file.exists(meta_file)) {
    stop("incorrect filepath for 'meta_file' provided")
  }
  if (is.null(out_dir)) {
    stop("'out_dir' is required")
  }
  if (!dir.exists(out_dir)) {
    stop("incorrect filepath for 'out_dir' provided")
  }
  
  ## load dataset intensity table
  ds_dat <- read.csv(ds_file, header = TRUE, stringsAsFactors = FALSE)
  ds_dat_cnames <- c("mz", "mzmin", "mzmax", "rt","rtmin", "rtmax","npeaks", "peakid", "pcs", "into")
  if (any(!ds_dat_cnames %in% colnames(ds_dat))) {
    stop("provided intensity table must contain columns: ", paste0(ds_dat_cnames, collapse = ", "))
  }
  
  ## load and order metadata by run order
  samples <- read.csv(meta_file, header = TRUE, stringsAsFactors = FALSE)
  samples <- samples[order(samples$run_order), ]
  
  ## for each pseudo chemical spectra, retain the intensity values from the sample with highest intensity for the corresponding peaks
  dat <- ds_dat[, -c(match(ds_dat_cnames, colnames(ds_dat)))]
  intens <- lapply(unique(ds_dat$pcs), FUN = getINTENSE, int_dat = ds_dat, dat = dat)
  intens <- do.call("rbind", intens)
  ds <- ds_dat[ , ds_dat_cnames[which(ds_dat_cnames != "into")]]
  ds <- cbind(ds, intens)
  
  object <- new(
    "massFlowAnno",
    filepath = ds_file,
    samples = samples,
    data = ds_dat,
    ds = ds
  )
  message("A 'massFlowAnno' object was succesfully built with ", nrow(samples), " samples.")
  return(object)
}

