#' Calculate Pearson correlations of a binned matrix
#' 
#' This function is used to generate a list object to be passed to getABSignal
#'
#' @param binmat      A binned matrix list object from getBinMatrix
#' 
#' @return    A list object to pass to getABSignal
#' 
#' @import    GenomicRanges
#' 
#' @export 

getCorMatrix <- function(binmat) {
  #Calculate correlations
  message("Calculating correlations...")
  binmat.cor <- cor(t(binmat$x))
  gr.cor  <- binmat$gr
  message("Done...")
  return(list(gr.cor=gr.cor, binmat.cor=binmat.cor))
}
