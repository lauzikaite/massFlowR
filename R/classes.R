#' Chemical reference database template.
#'
#' @slot filepath A \code{character} specifying path to database csv file.
#'
#' @description \code{massFlowDB} class encapsulates information provided in the chemical reference database.
#' Database contains all of the detected peaks for every chemical in the database. Peaks must be grouped into spectra using \code{pickPEAKS()} function prior to building of the \code{massFlowDB} class object.
#' The final database file is a single csv file, which must contain columns: "peakid", "mz", "rt", "into", "peakgr", "chemid", "dbid", "dbname".
#' Database will be used to initialise \code{massFlowTemplate} class object.
#' Database "peakid" and "peakgr" values will retain unmodified in the template, while every sample in the study will be aligned with them.
#'
#' @slot db A \code{data.frame} containing database template. Template must contain
#' @return
#' @export
#'
#' @examples
setClass("massFlowDB",
         slots = c(
           filepath = "character",
           db = "data.frame"
         ))


#' Sample alignment and annotation template.
#'
#' @slot filepath A \code{character} specifying path to a csv file with details on study sample names and their acquisition (run) order.
#' @slot samples A \code{data.frame} stores the "filepath" table.
#' @slot tmp A \code{data.frame} stores sample alignment and annotation template.
#'
#' @details \code{massFlowTemplate} object stores the sample alignment and annotation template.
#' Template is initiated using the first datafile in the study and (if provided) the chemical reference database (\code{massFlowDB} object).
#' To initiate the template, constructor function \code{\link{massFlowTemplate}} must be called, providing a "filepath" to the csv file with the details on the study samples. This csv file must contain columns "filepaths" and "run_order".
#' Each study sample must be processed with \code{\link{groupPEAKS}} function first in order to obtain spectral peak groups.
#' With every round of sample alignment & annotation, template is updated and an intermediate sample alignment output is written as a csv file.
#'
#' @return
#' @export
#'
#' @examples
setClass("massFlowTemplate",
         slots = c(
           filepath = "character",
           samples = "data.frame",
           tmp = "data.frame"))
