#' call mitochondrial variants from an MAlignments object 
#'
#' FIXME: figure out a way to reprocess extracted chrM/MT reads against rCRS,
#'        regardless of what reference they were originally aligned against.
#' 
#' @param mal         an MAlignments 
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
  if (!is(mal, "MAlignments")) stop("callMT needs an MAlignments to work.")
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
  res <- callVariants(qaVariants(tallyVariants(mal@bam, pars)), calling=filters)
  sampleNames(res) <- gsub(paste0(".", mtGenome), "", 
                           gsub("\\.bam", "", basename(mal@bam)))
  res$PASS <- apply(softFilterMatrix(res), 1, all) == 1
  res <- subset(res, totalDepth >= total.count)
  res$VAF <- altDepth(res) / totalDepth(res)
  genome(res) <- mtGenome
  mvr <- MVRanges(res, coverage(mal))
  if (rCRS == TRUE & (!mtGenome %in% c("GRCh38","hg38"))) {
    mvr <- rCRS(mvr)
  }
  names(mvr) <- mtHGVS(mvr)
  return(mvr)
}
