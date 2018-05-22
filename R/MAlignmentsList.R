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
  gal <- GenomicAlignments:::GAlignmentsList(...)
  metadata(gal)$coverage <- sapply(list(...), function(x) attr(x, "coverage"))
  if (is.null(names(gal))) {
    message("Warning: creating an MAlignmentsList without element names.")
  }
  mall <- new("MAlignmentsList", gal)
  return(mall)
}


#' number of aligned fragments per element
#'
#' (co-opting a generic from chromVAR to get reads per element)
#' 
#' @param object  an MAlignmentsList
#' 
#' @return        the fragments (length) for each element
#'
#' @import        chromVAR
#'
#' @export
setMethod("getFragmentsPerSample", signature(object="MAlignmentsList"),
          function(object) sapply(object, length))


#' estimated read coverage for each element (cached on load)
#' 
#' (co-opting a generic from `IRanges`)
#'
#' @param x   an MAlignmentsList
#' 
#' @return    estimated coverage (numeric vector)
#'
#' @import IRanges
#' 
#' @export
setMethod("coverage", signature(x="MAlignmentsList"),
          function(x) {
            covg <- metadata(x)$coverage
            # subset by sample name if sample names available 
            if (!is.null(names(x)) & !is.null(names(covg))) {
              covg <- covg[names(x)] 
            }
            return(covg)
          })


# resurrecting this from the `methods` package for use below...
setGeneric("getProperties", function(ClassDef) return(ClassDef))

#' element properties for an MAlignmentsList
#'
#' (co-opting a defunct `methods` function back into a generic)
#' 
#' @param ClassDef    an MAlignmentsList
#' 
#' @return            properties for each element, as a DataFrame
#'
#' @export
setMethod("getProperties", signature(ClassDef="MAlignmentsList"),
          function(ClassDef) {
            contig <- unique(seqlevels(ClassDef))
            basepairs <- seqlengths(ClassDef)[contig]
            dat <- DataFrame(reads=getFragmentsPerSample(ClassDef),
                             contig=rep(contig, length(ClassDef)),
                             basepairs=rep(basepairs, length(ClassDef)),
                             coverage=paste0(round(coverage(ClassDef)), "x"))
            if (!is.null(names(ClassDef))) rownames(dat) <- names(ClassDef)
            return(dat)
          })


#' display alignment collections with element summaries
#'
#' @param objects   an MAlignmentsList
#' 
#' @export
setMethod("show", signature(object="MAlignmentsList"),
          function(object) {
            cat("MAlignmentsList object of length", length(object), "\n")
            cat("-------\n", sep = "")
            cat("Element summary by getProperties():\n")
            show(getProperties(object))
            cat("-------\n", sep = "")
            cat("seqinfo: ", summary(seqinfo(object)), "\n", sep = "")
          })
