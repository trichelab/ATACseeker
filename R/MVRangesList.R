#' like a VRangesList, but for mitochondria
#' 
#' @import VariantAnnotation
#' 
#' @exportClass MVRangesList
setClass("MVRangesList", contains="SimpleVRangesList")


#' wrap a VRangesList for mitochondrial use
#'
#' @param ...     the MVRanges elements forming the MVRangesList
#'
#' @return        the MVRangesList
#' 
#' @export
MVRangesList <- function(...) {
  new("MVRangesList", GenomicRangesList(...), elementType = "MVRanges")
}


#' estimated read coverage (recorded from the called MAlignments)
#' 
#' @param x   an MVRangesList
#' 
#' @return    estimated coverage from the called MAlignments[List]
#'
#' @export
setMethod("coverage", signature(x="MVRangesList"), 
          function(x) sapply(x, coverage))


#' simple helper to retrieve PASS'ing, coding variants from an MVRangesList
#'
#' @param x   an MVRangesList
#'
#' @return    subsets of MVRangesList elements passing filters in coding regions
#'
#' @import    Biostrings
#'
#' @export
setMethod("encoding", signature(x="MVRangesList"), 
          function(x) MVRangesList(lapply(x, encoding)))


#' display variant calls with overall read coverage estimate
#'
#' @param x   an MVRangesList
#'
#' @import S4Vectors
#' 
#' @export
setMethod("show", signature(object="MVRangesList"),
          function(object) {
            callNextMethod()
            coverages <- paste0(round(unname(sapply(object, coverage))), "x")
            cat(S4Vectors:::labeledLine("coverage", coverages), "\n")
          })
