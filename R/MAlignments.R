#' wraps a GAlignments with information about coverage and its target BAM file
#' 
#' @import GenomicRanges
#' 
#' @exportClass MAlignments
setClass("MAlignments",
         # for calculating coverage: 
         representation(readLength="integer", 
                        genomeSize="integer",
                        bam="character",
                        bai="character",
                        gr="GRanges"),
         contains="GAlignments")

#' wrap a GAlignments for easier stats
#'
#' @param gal         a GAlignments
#' @param gr          a GRanges (mtRanges)
#' @param bam         a bam filename
#' @param bai         a bai filename
#' @param readLength  a number (read length)
#' @param genomeSize  a number (number of bases in the host genome)
#'
#' @return an MAlignments 
#' 
#' @import GenomicAlignments
#' 
#' @export
MAlignments <- function(gal, gr, bam, bai, readLength=75, genomeSize=16571) {
  if (!is(gal, "GAlignments")) stop("gal must be a GAlignments. Exiting.")
  if (!is(gr, "GRanges")) stop("gr must be a GRanges. Exiting.")
  if (length(bam) > 1) stop("bam must be a string. Exiting.")
  if (length(bai) > 1) stop("bai must be a string. Exiting.")
  if (length(readLength) > 1) stop("readLength must be an integer. Exiting.")
  if (length(genomeSize) > 1) stop("genomeSize must be an integer. Exiting.")
  new("MAlignments", gal, 
      readLength=as.integer(readLength), 
      genomeSize=as.integer(genomeSize),
      bam=bam,
      bai=bai,
      gr=gr)
}

#' estimated read coverage
#' 
#' @param x   an MAlignments
#' 
#' @return    estimated coverage (numeric)
#'
#' @export
setMethod("coverage", signature(x="MAlignments"),
          function(x) (length(x) * x@readLength) / x@genomeSize)

#' display alignment records with overall coverage estimate
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

#' retrieve the rest of the stored information in an MAlignments object
#'
#' @param x   an MAlignments
#' 
#' @return    a list of metadata (gr, bam, bai, readLength, genomeSize)
#'
#' @export
setMethod("metadata", signature(x="MAlignments"),
          function(x) {
            mdat <- list(gr=x@gr, 
                         bam=x@bam, 
                         bai=x@bai, 
                         readLength=x@readLength,
                         genomeSize=x@genomeSize)
            return(mdat)
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
            with(metadata(subject), BamViews(bam, bai, bamRanges=gr))
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
