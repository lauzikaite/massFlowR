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
suppressWarnings(library(dplyr))
rda_dir <- "/Users/el1514/Documents/Projects/DBannotate/NPC-RPOS-LCE-SubDB"
db_dir <- "/Users/el1514/Documents/Projects/massFlowR/packageDEV/db"

## ----db_table--------------------------------------------------------------
# rda_dir <- "absolute_path_to_rda_files"
# db_dir <- "absolute_path_to_output_directory"
db_table <- getDB(dir = rda_dir, out_dir = db_dir)

## ----db_table_knitr, echo = FALSE------------------------------------------
dt <- head(db_table)
knitr::kable(x = dt,
             format = "html",
             row.names = F) %>% 
  kableExtra::kable_styling(full_width = TRUE, bootstrap_options = c("striped"), position = "left")

