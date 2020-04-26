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
suppressWarnings(library(xcms))
rda_dir <- "~/Documents/Projects/DBannotate/NPC-RPOS-LCE-SubDB"
db_dir <- "~/Documents/Projects/massFlowR/packageDEV/db"
out_directory <- getwd()
url <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/massFlowR.html"

## ----db_table, message = FALSE, echo = FALSE, results = 'hide'-------------
buildDB(rda_dir = rda_dir, out_dir = out_directory)
db_table <- read.csv(file.path(out_directory, "DBtemplate.csv"))

## ----db_table_visible, eval=FALSE------------------------------------------
#  # rda_dir <- "absolute_path_to_rda_files"
#  # out_dir <- "absolute_path_to_output_directory"
#  buildDB(rda_dir = rda_dir, out_dir = out_directory)
#  db_table <- read.csv(file.path(out_directory, "DBtemplate.csv"))

## ----db_table_knitr, echo = FALSE------------------------------------------
db_table <- db_table[which(db_table$chemid %in% c(1,2,3)),]
dt <- head(db_table)

kableExtra::kable_styling(
  knitr::kable(x = dt,
           format = "html",
           row.names = F),
  full_width = TRUE, bootstrap_options = c("striped"), position = "left")

