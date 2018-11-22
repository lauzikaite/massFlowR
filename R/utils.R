## Unsorted utility functions

# cleanPEAKS ------------------------------------------------------------------------------------------------------
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
  peak <- dt_unique[rn, ]
  peak <- dt[which(dt$mz == peak$mz &
                     dt$rt == peak$rt), ]
  ## arrange by peakid and return the most intense
  peak <- peak[order(peak$into, decreasing = T), ]
  peak <- peak[1, ]
  return(peak)
}


# rbindCLEAN ------------------------------------------------------------------------------------------------------
#' @title Bind a list of data frames
#'
#' @description Function binds a list of data frames into a single data frame and removes inherited row names.
#'
#' @param ... list of dataframes
#'
#' @return Function returns a data frame.
#'
rbindCLEAN <- function(...) {
  rbind(..., make.row.names = F)
}


# scaleEDGES ------------------------------------------------------------------------------------------------------
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


# getCORmat -------------------------------------------------------------------------------------------------------
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
  setNames(as.data.frame(t(utils::combn(ind, 2, simplify = T))), nm = c("from", "to"))
}


# buildGRAPH ------------------------------------------------------------------------------------------------------
#' @title Build a correlation network of peaks
#'
#' @description Function build a correlation network of peaks using \emph{igraph} package functionality.
#' Only connections above the selected correlation threshold (\code{cor_thr}) will be retained in the network. 
#' 
#' @details Function requires a data.frame with unique peak-peak pairs and correlation coefficients describing their relationship.
#' Such data frame can be built using \code{\link{getCORmat}} function.
#' This function is used to build correlation networks using coefficients that represent either EIC correlation of peaks in a single LC-MS sample, or inter-sample intensity correlation.
#' 
#' @param pkg_cor \code{data.frame} containing a row for each unique peak-peak pair. Columns 'from' and 'to' specify peak identifiers, column 'weight' represents correlation coefficient.
#' @param cor_thr \code{numeric} specifying a threshold for correlation coefficient below which peak pairs are omitted from the network.
#' @param plot \code{logical} whether to save a PNG of the built network.
#' @param title \code{character} specifying name for the plot, required if PLOT = TRUE.
#' @param out_dir \code{character} specifying absolute path to output directory, where a PNG would be written.
#' 
#' @return Function returns community membership for each peak.
#'
buildGRAPH <-
  function(pkg_cor,
           cor_thr,
           plot = TRUE,
           title = NULL,
           out_dir = NULL
           ) {
    
    g <- igraph::graph_from_data_frame(pkg_cor, directed = FALSE)
    
    ## delete edges between peak pairs with cor < thr
    g <- igraph::delete.edges(g,
                              igraph::E(g)[[which(igraph::E(g)$weight < cor_thr)]])
    
    ## detect community structure in the network where non-correlated pairs are omitted
    coms <- igraph::cluster_label_prop(g)
    mem <- igraph::membership(coms)
    
    if (plot == TRUE) {
      grid::grid.newpage()
      grDevices::pdf(
        width = 8,
        height = 8,
        file = paste0(out_dir, "/", title, "_network.pdf")
      )
      ## scale edge thickness
      igraph::E(g)$weight <- scaleEDGES(igraph::E(g)$weight)
      
      ## make coordinates for vertices based on correlation: scale cor coef to maximise distance
      ## using Fruchterman-Reingold layout algorithm to prevent overlap
      coords <- igraph::layout_with_fr(g)
      igraph::plot.igraph(
        g,
        main = title,
        sub = paste("Cor threshold:", cor_thr),
        layout = coords
      )
      grDevices::dev.off()
    }
    
    # ## plot three mostly connected peaks
    # other_peaks <- V(g)$name[order(degree(g), decreasing = T)[-c(1:3)]]
    # top_peaks <- V(g)$name[order(degree(g), decreasing = T)[c(1:3)]]
    
    return(mem)
  }
