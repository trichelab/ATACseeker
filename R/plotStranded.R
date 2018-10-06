#' convenience function
#'
#' @param x     coverage
#' @param gr    a single GRange
#' @param ...   args to pass along to plotCoverage
#' 
#' @return    wrapper for plotting stranded coverage plots
#'
#' @import GenomicRanges
#'
#' @export
plotStranded <- function(x, gr, ...) plotCoverage(x, gr, ..., stranded=TRUE)

