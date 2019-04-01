#' massFlowR: a package for pre-processing of high-resolution, untargeted, centroid LC-MS data
#'
#' massFlowR annotates and aligns structurally-related spectral peaks across LC-MS experiment samples.
#' Its pipeline consists of three stages:
#'
#' \itemize{
#' \item Individual samples processing
#' \item Multiple sample alignment
#' \item Alignment validation and peak filling
#' }
#'
#' @import foreach
#' @import doParallel
#' @importFrom grDevices dev.off pdf
#' @importFrom methods new show validObject
#' @importFrom stats cor median na.omit setNames
#' @importFrom utils packageVersion read.csv write.csv
"_PACKAGE"