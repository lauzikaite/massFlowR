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
url_p <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/processing.html"
url_a <- "https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/annotation.html"

## ---- out.width = "700px", echo = FALSE---------------------------------------
knitr::include_graphics("../man/figures/scheme.png")

