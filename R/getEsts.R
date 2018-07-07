#' Complexity estimation for ATACseq libraries.
#'
#' Estimate the complexity of a library or sample based on unique fragments
#' using Daley and Smith's implementation of Good-Toulmin rational function
#' approximation to solve the missing species problem.
#' 
#' Original functions in preseqR v2.0.1.1 for this were:
#' preseqR.rfa.curve and preseqR.rfa.species.accum.curve
#' 
#' The new functions as of the version 4.0.0 are:
#' 
#' ds.rSAC == preseqR.rfa.curve
#' ds.rSAC.bootstrap == preseqR.rfa.species.accum.curve
#' 
#' The new functions return generators that can have data passed to them
#' instead of returning a data frame as in version 2.0.1.1.
#'
#' @param xx      The fragments or sample of fragments
#' @param withCI  Have preseq compute 95 percent confidence intervals for plots?
#' @param ...     Other arguments to pass on to preseqR 
#' 
#' @return        A data frame with results
#' 
#' @import preseqR
#' 
#' @export 
getEsts <- function(xx, withCI=FALSE, ...) {  
  if (is(xx, "GRanges") | is(xx, "GAlignmentPairs")) xx <- getComplexity(xx)
  message("Estimating complexity (this can take a little while)...")
  if (!withCI) {
    res.fun <- suppressWarnings(ds.rSAC(xx, mt = 100, ...))
    message("Done generating complexity estimator...")
    message("Building library complexity curves...")
    res <- data.frame(reads = as.double(xx[, 1] %*% xx[, 2]) * 1:100,
                      frags = sapply(1:100, res.fun$f))
  } else {
    res.fun <- suppressWarnings(ds.rSAC.bootstrap(xx, mt = 100, times = 99, conf = 0.95, ...))
    message("Done generating complexity estimator...")
    message("Building library complexity curves...")
    res.ci.lb <- sapply(1:100, res.fun$lb)
    res.ci.ub <- sapply(1:100, res.fun$ub)
    res <- data.frame(reads = as.double(xx[, 1] %*% xx[, 2]) * 1:100,
                      frags = sapply(1:100, res.fun$f),
                      lb = res.ci.lb,
                      ub = res.ci.ub)
  }
  return(res)
}
