#' Plots A/B compartment estimates on a per chromosome basis
#'
#' Plot A/B compartments bins
#'
#' @param x      The GRanges object returned from getABSignal
#' @param main    Title for the plot
#' @param ylim      Y-axis limits (default is -1 to 1)
#' @param unitarize   Should the data be unitarized?
#' @param reverse    Reverse the sign of the PC values?
#' @param top.col    Top (pos. PC values) chromatin color to be plotted
#' @param bot.col    Bottom (neg. PC values) chromatin color to be plotted
#' 
#' @return    invisibly, the compartment estimates from the plot
#' 
#' @import    GenomicRanges
#' 
#' @export 

plotAB <- function(x, main="",ylim=c(-1, 1), unitarize=FALSE, reverse=FALSE,
                      top.col = "deeppink4", bot.col = "grey50"){
  if (unitarize){
    x <- .unitarize(x)
  }
  x <- as.numeric(x)
  if (reverse){
    x <- -x
  }
  
  n <- length(x)
  col <- rep(top.col, n)
  col[x<0] <- bot.col
  barplot(x, ylim=ylim, 
          bty="n", xlab="", ylab="",border=col, col=col, main=main)
}

#Unitarize function taken from minfi::compartments.R
# Author: Jean-Philippe Fortin
# May 6th 2015
.unitarize <- function(x, medianCenter = TRUE) {
  if (medianCenter) x <- x - median(x, na.rm = TRUE)
  bad <- is.na(x)
  x[!bad] <- x[!bad] / sqrt(sum(x[!bad]^2))
  n.bad <- sum(bad)
  if (n.bad > 0) {
    message(
      sprintf("[.unitarize] %i missing values were ignored.\n", n.bad))
  }
  x
}
