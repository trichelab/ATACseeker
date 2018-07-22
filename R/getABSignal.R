#' Calculate Pearson correlations of smoothed eigenvectors
#' 
#' This function is used to generate a list x to be passed to getABSignal
#'
#' @param x      A list object from getCorMatrix
#' @param k    Bin size for smoothing (default = 2)
#' @param iter    Number of iterations for moving average smoothing (default = 2)
#' 
#' @return    A list x to pass to getABSignal
#' 
#' @import    GRanges
#' @import    mixOmics
#' 
#' @export 

getABSignal <- function(x, k = 2, iter = 2){
  message("Calculating eigenvectors...")
  pc <- .getFirstPC(x$binmat.cor)
  message(paste0("Smoothing with a bin size of ", k, " and ", iter, " iterations..."))
  pc <- .meanSmoother(pc, k=k, iter=iter)
  message("Done smoothing...")
  return(list(pc=pc, gr=x$gr))
}

#Internal function
#Code modified from https://github.com/Jfortin1/scATAC_Compartments
.getFirstPC <- function (matrix, ncomp = 1) {  
  matrix <- t(scale(t(matrix), center = TRUE, scale = FALSE))
  if (ncomp > 1) {
    p.mat <- nipals(matrix, ncomp = ncomp)$p
  }
  else {
    p.mat <- nipals(matrix, ncomp = ncomp)$p[, 1]
    csums <- colSums(matrix, na.rm=TRUE)
    if (cor(csums, p.mat) < 0){
      p.mat <- -p.mat
    }
  }
  p.mat <- p.mat * sqrt(length(p.mat)) #Chromosome length normalization
  return(p.mat)
}

#Internal function
#Code modified from minfi
# Author: Jean-Philippe Fortin
# May 6th 2015

.meanSmoother <- function(x, k = 2L, iter = 2L, na.rm = TRUE) {
    n <- length(x)
    y <- rep(NA_real_, n)
    
    window.mean <- function(x, j, k, na.rm = na.rm){
      if (k >= 1) {
        return(mean(x[seq(j - (k + 1L), j + k)], na.rm = na.rm))
      } else {
        x[j]
      }
    }
    
    for (i in (seq(k + 1L, n - k))) {
      y[i] <- window.mean(x, i, k, na.rm)
    }
    for (i in seq_len(k)) {
      y[i] <- window.mean(x, i, i - 1L, na.rm)
    }
    for (i in seq(n - k + 1L, n)) {
      y[i] <- window.mean(x, i, n - i, na.rm)
    }
    y
}
