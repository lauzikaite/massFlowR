## Unsorted utility functions.

#' @title Clean peak table from duplicating peaks
#' 
#' @description Function for every duplicated peak returns a single entry. 
#' Peak duplication can appear because of:
#' \itemize{
#' \item artefacts in centWave peak table (identical mz/rt) (inside \code{\link{pickPEAKS}} function)
#' \item different datafile-of-interest peaks matching the same template peaks (identical mz/rt/peakid) (inside \code{\link{addDOI}} function)
#' }
#' 
#' @param rn \code{numeric} specifying a single row number (rowname) of the peak table.
#' @param dt_unique \code{data.frame} with unique peaks' \emph{m/z} and \emph{rt} values, with as many rows, as there are unique peaks.
#' @param dt \code{data.frame} object containing the peak table from which duplicated peaks should be removed.
#'
#' @return Function returns a single peak entry for the unique \emph{m/z} and \emph{rt} combination.
#'
cleanPEAKS <- function(rn, dt_unique, dt) {
  ## extract full peak table for the corresponding peak
  peak <- dt_unique[rn,]
  peak <- dt[which(dt$mz == peak$mz &
                     dt$rt == peak$rt),]
  ## arrange by peakid and return the most intense
  peak <- peak[order(peak$into, decreasing = T),]
  peak <- peak[1,]
  return(peak)
}

#' @title Bind a list of data frames
#' 
#' @description Function binds a list of data frames into a single data frame and removes inherited row names.
#'
#' @param ... 
#'
#' @return Function returns a data frame.
#' 
rbindCLEAN <- function(...) {
  rbind(..., make.row.names = F)
}

#' @title Scale correlation coefficients
#' 
#' @description Function scales correlation coefficients to make a clear Fruchterman-Reingold graph.
#'
#' @param x \code{numeric} specifying correlation coefficients of a graph.
#' @param from \code{numeric} specifying the lowest value to which graph weights should be scaled to.
#' @param to \code{numeric} specifying the highest value to which graph weights should be scaled to.
#'
#' @return Functions returns scaled graph's weights in \code{numeric} format. 
#' 
scaleEDGES <- function(x, from = 0.01, to = 10) {
  (x - min(x)) / max(x - min(x)) * (to - from) + from
}

#' @title Build a correlation matrix between peaks-of-interest
#' 
#' @description Function builds an matrix for pairs between provided peaks-of-interest.
#' Each row is a unique peak-pair with columns 'from' and 'to' specifying peak indeces.
#' 
#' @param ind \code{numeric} with indeces of peaks-of-interest.
#'
#' @return Function returns a \code{data.frame} with unique peak-pairs.
#' 
getCORmat <- function(ind) {
  setNames(as.data.frame(t(combn(ind, 2, simplify = T))), nm = c("from", "to")) 
}

