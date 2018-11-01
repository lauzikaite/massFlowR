## ----setup, echo = FALSE, results = "asis"---------------------------------
BiocStyle::markdown() 
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
Biocpkg <- function (pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}
options(kpb.suppress_noninteractive = TRUE) ## supressing progress bar

## ----libraries, message = FALSE, echo = FALSE------------------------------
## load libraries quietly to avoid printing messages in the vignette
suppressWarnings(library(massFlowR))
suppressWarnings(library(xcms))
suppressWarnings(library(dplyr))

## ----installation, eval = FALSE--------------------------------------------
#  if (!require("devtools")) install.packages("devtools")
#  devtools::install_github("lauzikaite/massflowR")

## ----data_import-----------------------------------------------------------
library(massFlowR)
## Get the full path to the CDF files
# faahKO_files <- system.file(c('cdf/WT/wt15.CDF', 'cdf/WT/wt16.CDF'), package = "faahKO")
faahKO_files <- dir(system.file("cdf/WT", package = "faahKO"), full.names = T)
faahKO_files

## ----params----------------------------------------------------------------
## xcms parameters for peak-picking
cwt_param <- xcms::CentWaveParam(ppm = 25,
                                 snthresh = 10,
                                 noise = 1000,
                                 prefilter =  c(3, 100),
                                 peakwidth = c(30, 80))

