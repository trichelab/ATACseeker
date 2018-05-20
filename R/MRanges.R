#' like a VRanges, but for mitochondria
#' 
#' @import VariantAnnotation
#' 
#' @exportClass MVRanges
setClass("MVRanges", 
         representation(coverage="numeric"),
         contains="VRanges")

#' wrap a VRanges for mitochondrial use
#'
#' @param   vr    the VRanges
#' @param   covg  estimated coverage
#'
#' @return        an MVRanges
#' 
#' @export
MVRanges <- function(vr, coverage) new("MVRanges", vr, coverage=coverage)

#' estimated read coverage (recorded from the called MAlignments)
#' 
#' @param x   an MVRanges
#' 
#' @return    estimated coverage (numeric) from the called MAlignments
#'
#' @export
setMethod("coverage", signature(x="MVRanges"), function(x) x@coverage)

#' display variant calls with overall read coverage estimate
#'
#' @param x   an MVRanges
#' 
#' @export
setMethod("show", signature(object="MVRanges"),
          function(object) {
            callNextMethod()
            cat(paste0("  coverage: ~", round(coverage(object)), "x"), "\n")
          })
