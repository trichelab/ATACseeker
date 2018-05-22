#' call mitochondrial variants from an MAlignments object 
#'
#' FIXME: figure out a way to reprocess extracted chrM/MT reads against rCRS,
#'        regardless of what reference they were originally aligned against.
#' 
#' @param mal         an MAlignments (or, potentially, an MAlignmentsList) 
#' @param p.lower     lower bound on binomial probability for a variant (0.1)
#' @param read.count  minimum alt read depth required to support a variant (2)
#' @param total.count minimum total read depth required to keep a variant (10)
#' @param rCRS        lift to rCRS if not already hg38/GRCh38? (FALSE) 
#'
#' @import gmapR
#' @import VariantTools
#' @import GenomicAlignments
#' @import GmapGenome.Hsapiens.hg19.chrM
#'
#' @export
callMT <- function(mal, p.lower=.1, read.count=2L, total.count=10L, rCRS=FALSE){

  if (!is(mal, "MAlignments") & !is(mal, "MAlignmentsList")) {
    stop("callMT needs an MAlignments or MAlignmentsList to call variants.")
  } else if (is(mal, "MAlignmentsList")) { 
    message("Variant-calling an MAlignmentsList. This may melt your machine.")
    return(VRangesList(lapply(mal, callMT, parallel=parallel)))
  }

  mtChr <- seqlevelsInUse(mal)
  mtGenome <- unique(genome(mal))
  gmapGenome <- paste("GmapGenome", "Hsapiens", mtGenome, mtChr, sep=".")
  requireNamespace(gmapGenome)
  try(attachNamespace(gmapGenome), silent=TRUE)
  genome(mal) <- paste(mtGenome, mtChr, sep=".")
  isCircular(mal) <- FALSE # for variant calling
  pars <- TallyVariantsParam(get(gmapGenome), 
                             minimum_mapq=20L,
                             high_base_quality=20L,
                             ignore_duplicates=TRUE, 
                             read_length=as.integer(median(qwidth(mal))-1), 
                             which=as(seqinfo(mal)[mtChr], "GRanges"),
                             indels=TRUE)
  filters <- VariantCallingFilters(read.count, p.lower)
  tallied <- tallyVariants(asBam(mal), pars)
  QAed <- qaVariants(tallied)
  res <- callVariants(QAed, calling.filters=filters)
  sampleNames(res) <- gsub(paste0(".", mtGenome), "", 
                           gsub("\\.bam", "", basename(asBam(mal))))
  res$PASS <- apply(softFilterMatrix(res), 1, all) == 1
  res <- subset(res, totalDepth(res) >= total.count)
  res$VAF <- altDepth(res) / totalDepth(res)
  genome(res) <- mtGenome
  mvr <- MVRanges(res, coverage(mal))
  if (rCRS == TRUE) mvr <- rCRS(mvr)
  names(mvr) <- mtHGVS(mvr)
  return(mvr)
}
