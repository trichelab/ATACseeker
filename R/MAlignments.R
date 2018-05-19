#' wraps a GAlignments with information about coverage and its target BAM file
#' 
#' @import GenomicRanges
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
#' @return an MAlignments 
#' 
#' @import GenomicAlignments
#' 
#' @export
MAlignments <- function(gal, bam) { 
  if (!is(gal, "GAlignments")) stop("gal must be a GAlignments. Exiting.")
  if (length(bam) > 1) stop("bam must be a string. Exiting.")
  new("MAlignments", gal, bam=bam)
}

#' estimated read coverage
#' 
#' @param x   an MAlignments
#' 
#' @return    estimated coverage (numeric)
#'
#' @export
setMethod("coverage", signature(x="MAlignments"),
          function(x) {
            readLength <- median(qwidth(x)) - 1
            genomeSize <- width(as(seqinfo(x)[seqlevelsInUse(x)], "GRanges"))
            return((length(x) * readLength) / genomeSize)
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
#' @import Rsamtools
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
#' @param x           an MAlignments
#' @param p.lower     lower bound on variant probability to call (0.1)
#' @param p.error     estimated error probability for a variant (0.001)
#' @param read.count  minimum number of reads required to support a variant (10)
#' 
#' @return            a VRanges object
#'
#' @import VariantTools
#'
#' @export
setMethod("callVariants", signature(x="MAlignments"), 
          function(x, 
                   p.lower=0.1, 
                   p.error=0.001, 
                   read.count=10) {
            # replace with direct calling ASAP 
            callMT(x, p.lower, p.error, read.count) 
          })
