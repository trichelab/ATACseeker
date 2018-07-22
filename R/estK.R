#' Estimate a value of k to use with RUV
#' 
#' This function is used to estimate a suitable value of k for use in RUV normalization
#'
#' @param dat      A data matrix where loci are rows and samples are columns
#' @param mod    A model matrix used to fit the data (variables of interest)
#' @param method    The method used to estimate k ("leek" (default - asymptotic) or "be" (permutation-based))
#' @param vfilter    Specify the number of genes to filter based on ranked variance across samples (default is to use all genes)
#' @param B    Specify the number of permutations to use if method = "be" (default = 20)
#' @param seed    Specify a random number for setting the seed if method = "be"
#' 
#' @return    An estimated value of k for using in RUVseq
#' 
#' @import    sva
#' 
#' @export 

estK <- function(dat, mod, method = "leek", vfilter = NULL, B = 20, seed = 777) {
  message("Estimating a value of k using surrogate variable analysis...")
  k <- sva::num.sv(dat = dat, mod = mod, method = method, vfilter = vfilter, B = B, seed = seed)
  if (k == 0) {
    message("k is estimated to be 0. You may not need to remove unwanted variance...")
    message("It is also possible that something failed and you will need to manually specify k for RUV...")
  }
  else {
    message(paste0("k is estimated to be: ", k))
  }
  return(k)
}