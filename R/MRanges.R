#' like a VRanges, but for mitochondria
#' 
#' @import VariantAnnotation
#' 
#' @exportClass MRanges
setClass("MRanges", 
         representation(coverage="numeric"),
         contains="VRanges")

#' wrap a VRanges for mitochondrial use
#'
#' @param   vr    the VRanges
#' @param   covg  estimated coverage
#'
#' @return        an MRanges
#' 
#' @export
MRanges <- function(vr, coverage) new("MRanges", vr, coverage=coverage)

#' estimated read coverage (recorded from the called MAlignments)
#' 
#' @param x   an MRanges
#' 
#' @return    estimated coverage (numeric) from the called MAlignments
#'
#' @export
setMethod("coverage", signature(x="MRanges"), function(x) x@coverage)

#' display variant calls with overall read coverage estimate
#'
#' @param x   an MRanges
#' 
#' @export
setMethod("show", signature(object="MRanges"),
          function(object) {
            callNextMethod()
            cat(paste0("  coverage: ~", round(coverage(object)), "x"), "\n")
          })
