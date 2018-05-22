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
  if (is.null(names(gal))) warning("MAlignmentsList with no element names!")
  mall <- new("MAlignmentsList", gal)
  metadata(mall)$Summary <- Summary(mall)
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
#' @import    IRanges
#' 
#' @export
setMethod("coverage", signature(x="MAlignmentsList"),
          function(x) Summary(x)[, "coverage"])


#' estimated read length for each element (cached on load)
#' 
#' (co-opting a generic from `S4Vectors`)
#'
#' @param x   an MAlignmentsList
#' 
#' @return    estimated coverage (numeric vector)
#'
#' @import    S4Vectors
#' 
#' @export
setMethod("runLength", signature(x="MAlignmentsList"),
          function(x) sapply(x, runLength))


#' summary of an MAlignmentsList
#'
#' (generic defined in base, no less)
#' 
#' @param x    an MAlignmentsList
#' 
#' @return     a DataFrame
#'
#' @export
setMethod("Summary", signature(x="MAlignmentsList"),
          function(x) {

            # use cached version, if available
            if ("Summary" %in% names(metadata(x))) { 
              dat <- metadata(x)$Summary
              if (!is.null(names(x))) dat <- dat[names(x),]
            } else { 
              dat <- DataFrame(
                reads=sapply(x, length),
                readLength=sapply(x, runLength),
                genomeSize=sapply(x, yieldSize)
              ) 
              dat$coverage <- with(dat, (reads * readLength) / genomeSize)
              if (!is.null(names(x))) rownames(dat) <- names(x)
            }
            
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
            cat("Summary(object):\n")
            show(Summary(object))
            cat("-------\n", sep = "")
            cat("seqinfo: ", summary(seqinfo(object)), "\n", sep = "")
          })
