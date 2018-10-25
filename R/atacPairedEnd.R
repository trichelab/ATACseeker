#' pull in (filter and assemble as an object) some or all of the "useful" reads
#' in a paired-end ATACseq run (i.e., the majority of ATACseq studies)
#'
#' @param bam       character string, the BAM file to parse
#' @param bamParams optional parameters to pass through to ScanBamParam (none)
#' @param which     optional GRanges with specific regions to extract (all)
#'
#' @return  a GenomicAlignmentPairs object
#'
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import Homo.sapiens
#' @import Mus.musculus
#' @import Rsamtools
#'
#' @examples
#' \dontrun{
#'
#' data(cytobands_hg19)
#' library(Homo.sapiens)
#' del7q <- byBand(cytobands_hg19)$chr7q22.1
#' galp <- atacPairedEnd('chr7.q10.CD49f_r1.hg19.bam', which=del7q)
#'
#' library(Mus.musculus)
#' chr6 <- GRanges('chr6', IRanges(1, seqlengths(Mus.musculus)['chr6']), '*')
#' galp <- atacPairedEnd("A2.mm10.unique.bam",
#'                       bamParams=properPairedEndAtacFilters(which=chr6))
#'
#' }
#'
#' @export
atacPairedEnd <- function(bam, bamParams=NULL, which=NULL, ...) {


  if(is.null(bamParams)) {
    if (is.null(which)) {
      bf <- BamFile(bam) # pull in the index too 
      sl <- grep("(_|chrM|MT)", invert=TRUE, value=TRUE, seqlevels(bf))
      which <- as(seqinfo(bf)[sl], "GRanges")
    } 
    bamParams <- properPairedEndAtacFilters(which=which, ...)
  }

  readGAlignmentPairs(bam, param=bamParams)

}

## not exported; used for preseq estimation
properPairedEndAtacFilters <- function(which, ...) {

  ScanBamParam(what=c("rname","strand","pos","isize","mapq"),
               flag=scanBamFlag(isProperPair=TRUE,
                                isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

## not exported; used for BAM filtering, also needs which(!chrM) filter (duh?)
uniquePairedEndAtacFilters <- function(which, ...) {

  ScanBamParam(what=c("rname","strand","pos","isize","mapq"),
               flag=scanBamFlag(isProperPair=TRUE,
                                isDuplicate=FALSE,
                                isNotPrimaryRead=FALSE,
                                isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

## for recovering spike-ins
spikeInFilters <- function(...) {

  ## grab the unmapped sequences for brute-force alignment to phiX spike-ins
  ScanBamParam(what=c("seq"), 
               flag=scanBamFlag(isUnmappedQuery=TRUE, 
                                isNotPassingQualityControls=FALSE, 
                                ...))
}
