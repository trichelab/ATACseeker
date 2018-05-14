#' wraps a GAlignments with coverage information
#' 
#' @exportClass MAlignments
setClass("MAlignments",
         # for calculating coverage: 
         representation(readLength="integer", 
                        genomeSize="integer"),
         contains="GAlignments")

#' wrap a GAlignments for easier stats
#'
#' @param gal         a GAlignments
#' @param readLength  a number (read length)
#' @param genomeSize  a number (number of bases in the mitochondrial genome)
#'
#' @return an MAlignments 
#' 
#' @import GenomicAlignments
#' 
#' @export
MAlignments <- function(gal, readLength=75, genomeSize=16571) {
  if (length(readLength) > 1) stop("readLength must be an integer. Exiting.")
  if (length(genomeSize) > 1) stop("genomeSize must be an integer. Exiting.")
  new("MAlignments", gal, 
      readLength=as.integer(readLength), 
      genomeSize=as.integer(genomeSize))
}

setMethod("coverage", signature(x="MAlignments"),
          function(x) (length(x) * x@readLength) / x@genomeSize)

setMethod("show", signature(object="MAlignments"),
          function(object) {
            callNextMethod()
            cat("  -------\n")
            cat(paste0("  ", round(coverage(object)), 
                       "x approximate mtDNA read coverage."), "\n")
          })
