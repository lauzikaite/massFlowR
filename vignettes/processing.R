## ----setup, echo = FALSE, results = "asis"------------------------------------
BiocStyle::markdown() 
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
Biocpkg <- function (pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}
## supressing progress bar
options(kpb.suppress_noninteractive = TRUE) 

## ----libraries, message = FALSE, echo = FALSE---------------------------------
## load libraries quietly to avoid printing messages in the vignette
suppressWarnings(library(massFlowR))
out_directory <- getwd()
url_m <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/massFlowR.html"
url_a <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/annotation.html"

## ----data_import--------------------------------------------------------------
library(massFlowR)
## Get the full path to the CDF files
faahKO_files <- dir(system.file("cdf/WT", package = "faahKO"), full.names = TRUE)

## ----data_import_knitr, echo = FALSE------------------------------------------
kableExtra::kable_styling(
  knitr::kable(x = faahKO_files,
             format = "html",
             row.names = FALSE, col.names = ""),
  full_width = TRUE, bootstrap_options = c("striped", "hover", "responsive"), position = "left")

## ----params-------------------------------------------------------------------
## xcms parameters for peak-picking
cwt_param <- xcms::CentWaveParam(ppm = 25,
                                 snthresh = 10,
                                 noise = 1000,
                                 prefilter =  c(3, 100),
                                 peakwidth = c(30, 80),
                                 mzdiff = 0)

## ----metadata-----------------------------------------------------------------
## define where processed datafiles should be written
# out_directory <- "absolute_path_to_output_directory"

## create metadata table with required columns 'filename', 'run_order' and 'raw_filepath'
metadata <-
  data.frame(
  filename = gsub(".CDF", "", basename(faahKO_files)),
  run_order = seq(length(faahKO_files)),
  raw_filepath = faahKO_files,
  stringsAsFactors = FALSE
  )
write.csv(metadata, file = file.path(out_directory, "metadata.csv"), row.names = FALSE)

## ----experiment_table, echo = FALSE-------------------------------------------
kableExtra::kable_styling(
  knitr::kable(x = metadata,
             format = "html",
             row.names = FALSE),
  full_width = TRUE, bootstrap_options = c("striped", "hover", "responsive"), position = "left")

## ----peak_grouping------------------------------------------------------------
## run peak detection and grouping for the listed faahKO files with two workers
groupPEAKS(file = file.path(out_directory, "metadata.csv"), out_dir = out_directory, cwt = cwt_param, ncores = 2)

## ----metadata_update----------------------------------------------------------
##  update previous metadata table and add paths to generated csv files
processed_files <- list.files(out_directory, "peakgrs.csv", full.names = TRUE)
metadata$proc_filepath <- processed_files
write.csv(metadata, file.path(out_directory, "metadata_grouped.csv"), row.names = FALSE)

## ----metadata_update_table, echo = FALSE--------------------------------------
kableExtra::kable_styling(
  knitr::kable(x = metadata,
             format = "html",
             row.names = FALSE),
  full_width = TRUE, bootstrap_options = c("striped", "hover", "responsive"), position = "left")

## ----buildTMP-----------------------------------------------------------------
## initiate template
template <- buildTMP(file = file.path(out_directory, "metadata_grouped.csv"), out_dir = out_directory, mz_err = 0.01, rt_err = 10, cutoff = 0.5)

## ----samples_slot, eval = FALSE-----------------------------------------------
#  ## review samples in the experiment using slot "samples"
#  template@samples

## ----samples_slot_knitr, echo = FALSE-----------------------------------------
## review alignment results for an individual sample, e.g. the second, using slot "data"
dt <- template@samples
kableExtra::kable_styling(
  knitr::kable(x = dt,
             format = "html",
             row.names = F),
  full_width = TRUE, bootstrap_options = c("striped", "hover", "responsive"), position = "left")

## ----tmp_slot, eval = FALSE---------------------------------------------------
#  ## review gnerated template using slot "tmp"
#  head(template@tmp, 20)

## ----tmp_slot_knitr, echo = FALSE---------------------------------------------
## review alignment results for an individual sample, e.g. the second, using slot "data"
dt <- head(template@tmp, 20)
kableExtra::kable_styling(
  knitr::kable(x = dt,
             format = "html",
             row.names = F),
  full_width = TRUE, bootstrap_options = c("striped", "hover", "responsive"), position = "left")

## ----alignPEAKS---------------------------------------------------------------
## align peaks across all remaining samples
template <- alignPEAKS(template, out_dir = out_directory, ncores = 2)

## ----data_slot, eval = FALSE--------------------------------------------------
#  ## review alignment results for an individual sample, e.g. the second, using slot "data"
#  head(template@data[[2]])

## ----data_slot_knitr, echo = FALSE--------------------------------------------
## review alignment results for an individual sample, e.g. the second, using slot "data"
dt <- head(template@data[[2]])
kableExtra::kable_styling(
  knitr::kable(x = dt,
             format = "html",
             row.names = F),
  full_width = TRUE, bootstrap_options = c("striped", "hover", "responsive"), position = "left")

## ----loadALIGNED1-------------------------------------------------------------
## get the absolute paths to the updated metadata file and the final template writen by "alignPEAKS" 
m_file <- file.path(out_directory, "aligned.csv")
tmp_file <- file.path(out_directory, "template.csv")

## load already aligned samples into a massFlowTemplate object
template <- loadALIGNED(file = m_file, template = tmp_file, mz_err = 0.01, rt_err = 10, cutoff = 0.5)

## ----loadALIGNED2-------------------------------------------------------------
## if some of the samples in the experiment have not been aligned yet, re-start alignment
template@samples$aligned

## ----loadALIGNED3-------------------------------------------------------------
if (any(template@samples$aligned == FALSE)) {
  template <- alignPEAKS(template, out_dir = out_directory, ncores = 2)
}

## ----realtime, eval = FALSE---------------------------------------------------
#  template <- buildTMP(file = file.path(out_directory, "metadata_grouped.csv"), out_dir = out_directory, mz_err = 0.01, rt_err = 10, cutoff = 0.5, realtime = TRUE)
#  template <- alignPEAKS(template, out_dir = out_directory, ncores = 2)

## ----loadALIGNED--------------------------------------------------------------
## get the absolute paths to the updated metadata file and the final template writen by "alignPEAKS" 
m_file <- file.path(out_directory, "aligned.csv")
tmp_file <- file.path(out_directory, "template.csv")

## initiate validation by first loading aligned samples into a massFlowTemplate object
template <- loadALIGNED(file = m_file, template = tmp_file)

## ----validPEAKS---------------------------------------------------------------
## Start validation using a massFlowTemplate object
template <- validPEAKS(template, out_dir = out_directory, ncores = 2, cor_thr = 0.5)

## ----fillPEAKS----------------------------------------------------------------
## Fill peaks using validated massFlowTemplate object
template <- fillPEAKS(template, out_dir = out_directory, ncores = 2)

## ----load_filled--------------------------------------------------------------
final_dt <- read.csv(file.path(out_directory, "filled_intensity_data.csv"), header = TRUE, stringsAsFactors = FALSE)

## ----ion_map, echo = FALSE----------------------------------------------------
ggplot2::ggplot(final_dt,
             ggplot2::aes(x = rt, y = mz, color = as.factor(pcs))) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
  ggplot2::xlab("Retention time") +
  ggplot2::ylab("m/z") +
  ggplot2::scale_color_viridis_d(name = "Pseudo chemical spectra") +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(size = 0.1)
  )

## ----teardown, message = FALSE, echo = FALSE----------------------------------
unlink(list.files(path = out_directory, pattern = ".csv", full.names = T))
unlink(list.files(path = out_directory, pattern = ".RDS", full.names = T))

