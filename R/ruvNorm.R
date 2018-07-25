#' Remove unwanted variance (RUV)
#' 
#' This function is a wrapper for RUVs from RUVSeq to remove unwanted variance using control genes.
#' Additionally, if the number of factors to remove is unknown, the value of k is estimated using
#' surrogate variable analysis. However, if the number of factors is known a priori, then k can be
#' specified.
#'
#' @param x      A data matrix where loci are rows and samples are columns
#' @param mod    A model matrix of variables to fit the data
#' @param cIdx    A vector of loci/spike-ins in the data matrix to be used as negative controls (default = all loci)
#' @param k    The number of factors of unwanted variation to be estimated/removed from the data (default = estimate)
#' @param scIdx    A matrix of sample replicates to be used to estimate/remove unwanted variation
#' @param round    Measures are rounded to form pseudo-counts
#' @param epsilon    An offset to be added to avoid taking log of 0 (default = 1)
#' @param tolerance    Selection tolerance for singular values to be considered positive (default = 1e-8)
#' @param isLog    Whether the input matrix is already log-transformed (default = FALSE)
#' @param ...    Other arguments to be passed to estK
#' 
#' @return    A list object of the samples-by-factors matrix of estimated unwanted variation and the normalized counts (removed unwanted variation)
#' 
#' @import    RUVSeq
#' 
#' @export 

ruvNorm <- function(x, mod, cIdx, k = "estimate", scIdx, round = TRUE, epsilon = 1, tolerance = 1e-8, isLog = FALSE, ...) {
  if (k == "estimate") {
    if (!class(mod) == "matrix") {
      stop("mod needs to be a model.matrix object...")
    }
    k <- estK(dat = x, mod = mod, ...)
  }
  #TODO: make the scIdx generation automated
  message("Removing unwanted variance...")
  ruv.out <- RUVs(x = x, cIdx = rownames(x), k = k, scIdx = scIdx, round = round, epsilon = epsilon, tolerance = tolerance, isLog = isLog)
  return(ruv.out)
}
