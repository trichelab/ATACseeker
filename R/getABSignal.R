#' Calculate Pearson correlations of smoothed eigenvectors
#' 
#' This function is used to generate a list x to be passed to getABSignal
#'
#' @param x      A list object from getCorMatrix
#' @param k    Value of k for smoothing (default = 2)
#' @param iter    Number of iterations for moving average smoothing (default = 2)
#' 
#' @return    A list x to pass to getABSignal
#' 
#' @import    GenomicRanges
#' @import    mixOmics
#' 
#' @export 

getABSignal <- function(x, k = 2, iter = 2){
  message("Calculating eigenvectors...")
  pc <- .getFirstPC(x$binmat.cor)
  message(paste0("Smoothing with a k of ", k, " for ", iter, " iterations..."))
  pc <- .meanSmoother(pc, k=k, iter=iter)
  message("Done smoothing...")
  gr <- x$gr
  gr$pc <- pc
  gr$compartments <- .extractOpenClosed(pc)
  return(gr)
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

.meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
  meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
    n <- length(x)
    y <- rep(NA,n)
    
    window.mean <- function(x, j, k, na.rm=na.rm){
      if (k>=1){
        return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
      } else {
        return(x[j])
      }    
    }
    
    for (i in (k+1):(n-k)){
      y[i] <- window.mean(x,i,k, na.rm)
    }
    for (i in 1:k){
      y[i] <- window.mean(x,i,i-1, na.rm)
    }
    for (i in (n-k+1):n){
      y[i] <- window.mean(x,i,n-i,na.rm)
    }
    y
  }
  
  for (i in 1:iter){
    x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
  }
  x
}

#Internal function
#Code modified from minfi
# Author: Jean-Philippe Fortin
# May 6th 2015
.extractOpenClosed <- function(pc, cutoff = 0){
  ifelse(pc < cutoff, "open", "closed")
}
