#' wraps a GAlignments with information about coverage and its target BAM file
#' 
#' @import GenomicAlignments
#' 
#' @exportClass MAlignments
setClass("MAlignments",
         representation(bam="character",
                        runLength="numeric"), 
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
  if (length(bam) > 1) stop("bam must be a string naming a BAM file. Exiting.")
  new("MAlignments", gal, bam=bam, runLength=(median(qwidth(gal))-1))
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
            unname( (length(x) * runLength(x)) / yieldSize(x) )
          })


#' median read length
#' 
#' @param x   an MAlignments
#' 
#' @return    read length
#'
#' @export
setMethod("runLength", signature(x="MAlignments"),
          function(x) {
            return(x@runLength)
          })


#' genome size
#' 
#' @param x   an MAlignments
#' 
#' @return    genome size 
#'
#' @export
setMethod("yieldSize", signature(object="MAlignments"),
          function(object) {
            unname(seqlengths(object)[seqlevelsInUse(object)])
          })


#' summary of object: reads, readLength, genomeSize, coverage
#' 
#' specifically, 
#' c(length(x), runLength(x), seqlengths(x)[seqlevelsInUse(x)], coverage(x))
#' 
#' @param x   an MAlignments
#' 
#' @return    summary of x
#'
#' @export
setMethod("Summary", signature(x="MAlignments"),
          function(x) {
            c(reads=length(x),
              readLength=runLength(x),
              genomeSize=yieldSize(x),
              coverage=coverage(x))
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


#' fetch the filename for an MAlignments
#'
#' @param file      an MAlignments
#' 
#' @return          a BAM filename (string)
#'
#' @import          Rsamtools
#'
#' @export
setMethod("asBam", signature(file="MAlignments"),
          function(file) {
            return(file@bam)
          })


#' fetch the header from an MAlignments' original BAM 
#'
#' @param bamRanges an MAlignments
#' 
#' @return          a BamViews object
#'
#' @import          Rsamtools
#'
#' @export
setMethod("scanBamHeader", signature(files="MAlignments"),
          function(files) {
            return(scanBamHeader(asBam(files)))
          })


#' recreate the BamViews used to load an MAlignments
#'
#' @param bamRanges an MAlignments
#' 
#' @return          a BamViews object
#'
#' @import          Rsamtools
#'
#' @export
setMethod("BamViews", signature(bamRanges="MAlignments"),
          function(bamRanges) {
            bam <- bamRanges@bam
            bai <- paste0(bam, ".bai")
            gr <- as(seqinfo(bamRanges)[seqlevelsInUse(bamRanges)], "GRanges")
            BamViews(bam=bam, bai=bai, bamRanges=gr)
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
