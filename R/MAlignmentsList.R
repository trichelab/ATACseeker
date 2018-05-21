#' wraps a GAlignmentsList (made up of MAlignments) for nicer viewing
#' 
#' @import GenomicAlignments
#' 
#' @exportClass MAlignmentsList
setClass("MAlignmentsList", contains="GAlignmentsList")


#' wrap a GAlignmentsList for viewing
#'
#' @param ...         MAlignments
#'
#' @return            an MAlignments 
#' 
#' @import            GenomicAlignments
#' 
#' @export
MAlignmentsList <- function(...) {
  mall <- new("MAlignments", GAlignmentsList(...))
  message("Calculating and caching read coverage...") 
  metadata(mall)$coverage <- sapply(mall, coverage)
  return(mall)
}


#' estimated read coverage for each element (cached on load)
#' 
#' @param x   an MAlignmentsList
#' 
#' @return    estimated coverage (numeric vector)
#'
#' @export
setMethod("coverage", signature(x="MAlignmentsList"),
          function(x) {
            metadata(x)$coverage[colnames(x)]
          })


#' display alignment collections with read coverage estimates
#'
#' @param objects   an MAlignmentsList
#' 
#' @export
setMethod("show", signature(object="MAlignmentsList"),
          function(object) {
            cat("MAlignmentsList object of length ", length(object), "\n")
            cat("Summary of elements:", "\n")
            dat <- DataFrame(reads=sapply(object, length),
                             coverage=coverage(object))
            rownames(dat) <- names(object)
            show(dat)
          })
