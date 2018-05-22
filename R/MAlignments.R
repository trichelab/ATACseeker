#' wraps a GAlignments with information about coverage and its target BAM file
#' 
#' @import GenomicAlignments
#' 
#' @exportClass MAlignments
setClass("MAlignments",
         representation(bam="character"), 
         contains="GAlignments")


#' wrap a GAlignments for easier stats
#'
#' @param gal         a GAlignments
#' @param bam         a bam filename
#'
#' @return            an MAlignments 
#' 
#' @import            GenomicAlignments
#' 
#' @export
MAlignments <- function(gal, bam) { 
  if (!is(gal, "GAlignments")) stop("gal must be a GAlignments. Exiting.")
  if (length(bam) > 1) stop("bam must be a string. Exiting.")
  x <- new("MAlignments", gal, bam=bam)
  metadata(x)$runLength <- runLength(x)
  genomeSize <- seqlengths(x)[seqlevelsInUse(x)]
  metadata(x)$coverage <- with(metadata(x), (length(x)*runLength) / genomeSize)
  metadata(x)$Summary <- Summary(x)
  return(x)
}


#' estimated read coverage (cached when possible)
#' 
#' @param x   an MAlignments
#' 
#' @return    estimated coverage (numeric)
#'
#' @export
setMethod("coverage", signature(x="MAlignments"),
          function(x) {
            if ("coverage" %in% names(metadata(x))) {
              return(metadata(x)$coverage)
            } else { 
              genomeSize <- seqlengths(x)[seqlevelsInUse(x)]
              return( (length(x)*runLength(x)) / genomeSize )
            }
          })


#' median read length (cached when possible)
#' 
#' @param x   an MAlignments
#' 
#' @return    estimated coverage (numeric)
#'
#' @export
setMethod("runLength", signature(x="MAlignments"),
          function(x) { 
            if ("readLength" %in% names(metadata(x))) {
              return(metadata(x)$readLength)
            } else { 
              return(median(qwidth(x)) - 1)
            }
          })


#' Summary: reads, readLength, genomeSize, coverage
#' 
#' specifically, 
#' c(length(x), runLength(x), seqlengths(x)[seqlevelsInUse(x)], coverage(x))
#' 
#' @param x   an MAlignments
#' 
#' @return    estimated coverage (numeric)
#'
#' @export
setMethod("Summary", signature(x="MAlignments"),
          function(x) {
            if ("Summary" %in% names(metadata(x))) {
              res <- metadata(x)$Summary
            } else { 
              res <- c(reads=length(x),
                       readLength=runLength(x),
                       genomeSize=seqlengths(x)[seqlevelsInUse(x)],
                       coverage=coverage(x))
            }
            return(res)
          })


#' display alignment records with read coverage estimate
#'
#' @param x   an MAlignments
#' 
#' @export
setMethod("show", signature(object="MAlignments"),
          function(object) {
            callNextMethod()
            cat("  -------\n")
            cat(paste0("  ", round(coverage(object)), 
                       "x approximate read coverage."), "\n")
          })


#' recreate the BamViews needed to call variants 
#'
#' @param subject   an MAlignments
#' 
#' @return          a BamViews object
#'
#' @import          Rsamtools
#'
#' @export
setMethod("Views", signature(subject="MAlignments"),
          function(subject) {
            bam <- subject@bam
            bai <- paste0(bam, ".bai")
            bamRanges <- as(seqinfo(subject)[seqlevelsInUse(subject)],"GRanges")
            BamViews(bam=bam, bai=bai, bamRanges=bamRanges)
          })


#' call variants (shockingly enough) 
#'
#' @param x   an MAlignments
#' 
#' @return    a VRanges object
#'
#' @import    VariantTools
#'
#' @export
setMethod("callVariants", signature(x="MAlignments"), function(x) callMT(x))
