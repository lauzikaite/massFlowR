## ----setup, echo = FALSE, results = "asis"---------------------------------
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

## ----libraries, message = FALSE, echo = FALSE------------------------------
## load libraries quietly to avoid printing messages in the vignette
suppressWarnings(library(massFlowR))
suppressWarnings(library(xcms))
out_directory <- getwd()
url_m <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/massFlowR.html"
url_a <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/annotation.html"

## ----installation, eval = FALSE--------------------------------------------
#  if (!require("devtools")) install.packages("devtools")
#  devtools::install_github("lauzikaite/massflowR")

## ----data_import-----------------------------------------------------------
library(massFlowR)
## Get the full path to the CDF files
# faahKO_files <- system.file(c('cdf/WT/wt15.CDF', 'cdf/WT/wt16.CDF'), package = "faahKO"))
faahKO_files <- dir(system.file("cdf/WT", package = "faahKO"), full.names = T)
faahKO_files

## ----params----------------------------------------------------------------
## xcms parameters for peak-picking
cwt_param <- xcms::CentWaveParam(ppm = 25,
                                 snthresh = 10,
                                 noise = 1000,
                                 prefilter =  c(3, 100),
                                 peakwidth = c(30, 80),
                                 mzdiff = 0)

## ----metadata--------------------------------------------------------------
## define where processed datafiles should be written
# out_directory <- "absolute_path_to_output_directory"

## create metadata table with required columns 'filename', 'run_order' and 'raw_filepath'
metadata <-
  data.frame(
  filename = gsub(".CDF", "", basename(faahKO_files)),
  run_order = 1:length(faahKO_files),
  raw_filepath = faahKO_files
  )
write.csv(metadata, file = file.path(out_directory, "metadata.csv"), row.names = FALSE)

## ----experiment_table, echo = FALSE----------------------------------------
kableExtra::kable_styling(
  knitr::kable(x = metadata,
             format = "html",
             row.names = F),
  full_width = TRUE, bootstrap_options = "striped", position = "left")

## ----peak_grouping---------------------------------------------------------
## run peak detection and grouping for the listed faahKO files with two workers
groupPEAKS(file = file.path(out_directory, "metadata.csv"), out_dir = out_directory, cwt = cwt_param, ncores = 2)

## ----metadata_update-------------------------------------------------------
##  update previous metadata table
processed_files <- list.files(out_directory, "peakgrs.csv", full.names = T)
metadata$proc_filepath <- processed_files
write.csv(metadata, file.path(out_directory, "metadata_grouped.csv"), row.names = F)

## ----metadata_update_table, echo = FALSE-----------------------------------
kableExtra::kable_styling(
  knitr::kable(x = metadata,
             format = "html",
             row.names = F),
  full_width = TRUE, bootstrap_options = "striped", position = "left")

## ----buildTMP--------------------------------------------------------------
## initiate template
template <- buildTMP(file = file.path(out_directory, "metadata_grouped.csv"), out_dir = out_directory, mz_err = 0.01, rt_err = 10)

## ----massFlowTemplate------------------------------------------------------
## review gnerated template using slot "tmp"
head(template@tmp)

## review samples in the experiment using slot "samples"
template@samples

## ----alignPEAKS------------------------------------------------------------
## align peaks across all remaining samples
template <- alignPEAKS(template, out_dir = out_directory, ncores = 2)

## ----data_slot, eval = FALSE-----------------------------------------------
#  ## review alignment results for an individual sample, e.g. the second, using slot "data"
#  head(template@data[[2]])

## ----data_slot_knitr, echo = FALSE-----------------------------------------
## review alignment results for an individual sample, e.g. the second, using slot "data"
dt <- head(template@data[[2]])
kableExtra::kable_styling(
  knitr::kable(x = dt,
             format = "html",
             row.names = F),
  full_width = TRUE, bootstrap_options = c("striped"), position = "left")

## ----loadALIGNED-----------------------------------------------------------
## get the absolute paths to the updated metadata file and the final template
m_file <- file.path(out_directory, "aligned.csv")
tmp_file <- file.path(out_directory, "template.csv")

## initiate validation by first loading aligned samples into a massFlowTemplate object
template <- loadALIGNED(file = m_file, template = tmp_file, rt_err = 10)

## ----validPEAKS------------------------------------------------------------
## Start validation using a massFlowTemplate object
template <- validPEAKS(template, out_dir = out_directory, ncores = 2, cor_thr = 0.5)

## ----fillPEAKS-------------------------------------------------------------
## Fill peaks using validated massFlowTemplate object
template <- fillPEAKS(template, out_dir = out_directory, ncores = 2)

## ----teardown, message = FALSE, echo = FALSE-------------------------------
unlink(list.files(path = out_directory, pattern = ".csv", full.names = T))
unlink(list.files(path = out_directory, pattern = ".RDS", full.names = T))

