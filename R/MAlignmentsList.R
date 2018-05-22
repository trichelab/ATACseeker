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
  if (is.null(names(gal))) warning("This MAlignmentsList has no element names!")
  mall <- new("MAlignmentsList", gal)
  message("Caching element summary statistics...") 
  metadata(mall)$Summary <- Summary(mall)
  message("...done. Caching BAM file listings...")
  metadata(mall)$asBam <- asBam(mall)
  message("...done.")
  return(mall)
}


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
          function(x) Summary(x)[, "genomeCoverage"])


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


#' list the BAM files and their aligned genomes for MAlignmentsList elements
#' 
#' (co-opting a generic from `Rsamtools`)
#'
#' @param file  an MAlignmentsList
#' 
#' @return      BAM file summary for the MAlignmentsList 
#'
#' @import      S4Vectors
#' 
#' @export
setMethod("asBam", signature(file="MAlignmentsList"),
          function(file) {
            if ("asBam" %in% names(metadata(file))) {
              BAMs <- metadata(file)$asBam
            } else {
              BAMs <- DataFrame(BAM=sapply(file, asBam),
                                genome=unname(sapply(file, genome)))
              if (!is.null(names(file))) rownames(BAMs) <- names(file)
            } 
            if (!is.null(names(file))) BAMs <- BAMs[names(file), ] 
            return(BAMs)
          })


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
            if ("Summary" %in% names(metadata(x))) {
              dat <- metadata(x)$Summary
            } else {
              dat <- DataFrame(reads=sapply(x, length),
                               readLength=sapply(x, runLength),
                               genomeSize=sapply(x, yieldSize))
              if (!is.null(names(x))) rownames(dat) <- names(x)
              dat$genomeCoverage <- round(with(dat,
                                               (reads*readLength) / genomeSize))
            }
            # try to make the most of cached summaries... 
            if (!is.null(names(x))) dat <- dat[names(x), ] 
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
            cat("Summary(object) returns a ")
            show(Summary(object))
            cat("-------\n", sep = "")
            cat("seqinfo: ", summary(seqinfo(object)), "\n", sep = "")
          })
