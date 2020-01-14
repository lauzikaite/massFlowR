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
suppressWarnings(library(massFlowR))
db_file <- "../inst/testdata/database.csv"
ds_file <- "../inst/testdata/filled_intensity_data.csv"
meta_file <- "../inst/testdata/metadata.csv"
out_directory <- getwd()
massFlowR_url <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/massFlowR.html"
processing_url <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/processing.html"

## ----db_table_visible, eval=FALSE------------------------------------------
#  # rda_dir <- "path to rda files"
#  # out_directory <- "path to output directory"
#  buildDB(rda_dir = rda_dir, out_dir = out_directory)
#  db_table <- read.csv(file.path(out_directory, "database.csv"))

## ----db_table_knitr, echo = FALSE------------------------------------------
db_table <- read.csv(db_file)
db_table <- db_table[which(db_table$chemid %in% c(1,2,3)),]
kableExtra::kable_styling(
  knitr::kable(x = db_table,
           format = "html",
           row.names = F),
  full_width = TRUE, bootstrap_options = c("striped"), position = "left")

## ----massFlowAnno----------------------------------------------------------
# meta_file <- "path to metadata csv file"
# ds_file <-  "path to filled intensity data csv file"
# out_directory <- "path to output directory"
anno <- buildANNO(ds_file = ds_file, meta_file = meta_file, out_dir = out_directory)

## ----annotateDS------------------------------------------------------------
anno <- annotateDS(object = anno, db_file = db_file, out_dir = out_directory, mz_err = 0.01, rt_err = 10, ncores = 2)

## ----teardown, message = FALSE, echo = FALSE-------------------------------
unlink(list.files(path = out_directory, pattern = ".csv", full.names = TRUE))
unlink(list.files(path = out_directory, pattern = ".RDS", full.names = TRUE))

