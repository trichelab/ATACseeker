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

#' MVRangesList methods (centralized).
#'
#' @name      MVRangesList-methods
NULL


#' @rdname    MVRangesList-methods
#' 
#' @param x   an MVRangesList
#' 
#' @return    estimated coverage from the called MAlignments[List]
#'
#' @export
setMethod("coverage", signature(x="MVRangesList"), 
          function(x) sapply(x, coverage))


#' @rdname    MVRangesList-methods
#' 
#' @param object  an MVRangesList
#' @param value   a RangedSummarizedExperiment with identical colnames
#' 
#' @export
setReplaceMethod("counts", 
                 signature(object="MVRangesList", 
                           value="RangedSummarizedExperiment"),
                 function(object, value) {
                   if (!identical(names(object), colnames(value))) {
                     stop("Error: colnames(value) doesn't match names(object)!")
                   } else if (!"counts" %in% names(assays(value))) {
                     stop("Error: value must have an assay named `counts`!")
                   } else { 
                     metadata(object)$counts <- value[, names(object)]
                     return(object)
                   }
                 })


#' @rdname    MVRangesList-methods
#' 
#' @param object  an MVRangesList
#' 
#' @export
setMethod("counts", signature(object="MVRangesList"), 
          function(object) metadata(object)$counts)


#' @rdname    MVRangesList-methods
#'
#' @param object        an MVRangesList with an RSE in metadata(object)$counts
#' @param annotations   RangedSummarizedExperiment with motif matches for object
#'
#' @import chromVAR
#' 
#' @export
setMethod("computeDeviations", 
          signature(object="MVRangesList", 
                    annotations="RangedSummarizedExperiment"),
          function(object, annotations) { 
            if (!"counts" %in% names(metadata(object))) {
              stop("Error: metadata(object)$counts is currently empty.")
            } else {
              computeDeviations(object=counts(object), annotations=annotations)
            }
          })


#' @rdname    MVRangesList-methods
#' 
#' @param object  an MVRangesList
#' 
#' @export
setMethod("deviations", signature(object="MVRangesList"), 
          function(object) metadata(object)$deviations)


setAs(from="MVRangesList", to="chromVARDeviations", 
      function(from) deviations(from))


#' @rdname    MVRangesList-methods
#'
#' @param x   an MVRangesList
#'
#' @import    Biostrings
#'
#' @export
setMethod("encoding", signature(x="MVRangesList"), 
          function(x) MVRangesList(lapply(x, encoding)))


#' @rdname    MVRangesList-methods
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
            cat(S4Vectors:::labeledLine("coverage", coverages))
            if ("counts" %in% names(metadata(object))) {
              peaks <- nrow(metadata(object)$counts)
              cat(ifelse("bias" %in% names(rowData(counts(object))),
                  "Bias-corrected ", "Raw "))
              cat("fragment counts at", peaks, "peaks are available from",
                  "counts(object).\n")
            }
          })
