#' like a VRangesList, but for mitochondria
#' 
#' @import VariantAnnotation
#' @import Biostrings
#' @import S4Vectors
#' @import chromVAR
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
#' `counts`     returns fragment counts, if any, as a SummarizedExperiment.
#' `counts<-`   adds or updates fragment counts from a SummarizedExperiment.
#' `coverage`   returns the estimated coverage for each element in the list.
#' `deviations` returns deviations computed and stored from computeDeviations.
#' `encoding`   returns a subset of mutations in coding regions of mtDNA genes.
#' 
#' @param x            an MVRangesList (for some methods)
#' @param object       an MVRangesList (for other methods)
#' @param annotations  a RangedSummarizedExperiment with motif matches
#' @param value        a RangedSummarizedExperiment with matching colnames
#'
#' @name  MVRangesList-methods
NULL


#' @rdname    MVRangesList-methods
#' @export
setMethod("coverage", signature(x="MVRangesList"), 
          function(x) sapply(x, coverage))


#' @rdname    MVRangesList-methods
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
                     columns <- names(object)
                     metadata(object)$counts <- filterPeaks(value[, columns])
                     return(object)
                   }
                 })


#' @rdname    MVRangesList-methods
#' @export
setMethod("counts", signature(object="MVRangesList"), 
          function(object) metadata(object)$counts)


#' @rdname    MVRangesList-methods
#' @export
setMethod("encoding", signature(x="MVRangesList"), 
          function(x) MVRangesList(lapply(x, encoding)))


#' @rdname    MVRangesList-methods
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
