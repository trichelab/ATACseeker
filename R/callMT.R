#' call mitochondrial variants from the results of ATACseeker::getMT() 
#'
#' call mitochondrial variants against the reference sequence they aligned to
#' the appropriate cutoffs will depend on the estimated copy number of mtDNA
#' in the cell type assayed; leukocytes usually have between 100-500, while
#' hepatocytes can have thousands. Therefore different settings make sense
#' for different cell types. The defaults call a mutation with ~ 2% VAF.
#'
#' FIXME: figure out a way to reprocess extracted chrM/MT reads against rCRS,
#'        regardless of what reference they were originally aligned against.
#' 
#' @param mtReads     mitochondrial reads, with bamViews, or a BAM filename
#' @param p.lower     lower bound on binomial probability for a variant (0.1)
#' @param p.error     error probability (influences the minimum VAF; 0.001)
#' @param read.count  minimum read depth required to support a variant (2)
#' @param ...         any other arguments to pass (currently ignored)
#'
#' @import gmapR
#' @import VariantTools
#' @import GmapGenome.Hsapiens.hg19.chrM
#' @import GmapGenome.Hsapiens.hg38.chrM
#' @import GmapGenome.Hsapiens.GRCh38.MT
#'
#' @export
callMT <- function(mtReads, p.lower=0.1, p.error=0.001, read.count=2L, ...) { 
  
  if (is(mtReads, "character")) mtReads <- getMT(mtReads)
  if (!is(mtReads, "MAlignments")) {
    stop("Need an 'enhanced' MAlignments result from getMT() for this to work.")
  }
  mtGenome <- unique(genome(mtReads))
  mtChr <- unique(names(genome(mtReads)))
  gmapGenome <- paste("GmapGenome", "Hsapiens", mtGenome, mtChr, sep=".")
  requireNamespace(gmapGenome)
  try(attachNamespace(gmapGenome), silent=TRUE)
  bam <- bamPaths(mtView(mtReads))
  readLen <- mtReads@readLength
  whichRanges <- bamRanges(BamViews(bam, bai, bamRanges = mtRange))

  # FIXME: use a mask= argument to black out hypervariable region?
  tally.param <- TallyVariantsParam(get(gmapGenome),
                                    minimum_mapq=20,
                                    high_base_quality=20L,
                                    read_length=readLen, 
                                    ignore_duplicates=TRUE, 
                                    which=whichRanges,
                                    indels=TRUE)
  tallies <- tallyVariants(bam, tally.param) 
  qa.variants <- qaVariants(tallies)
  calling.param <- VariantCallingFilters(read.count=read.count,
                                         p.lower=p.lower,
                                         p.error=p.error)
  res <- callVariants(qa.variants, calling.filters=calling.param)
  sampleNames(res) <- gsub(paste0(".", mtGenome), "", 
                           gsub("\\.bam", "", 
                                basename(bamPaths(attr(mtReads, "mtView")))))
  res$PASS <- apply(softFilterMatrix(res), 1, all) == 1
  res <- res[rev(order(res$PASS, totalDepth(res)))]
  res$VAF <- altDepth(res) / totalDepth(res)
  return(res)
}

# hg19, hg38, and GRCh38 mitochondrial genome creation for gmapR
# (any genome with a FASTA of the MT contig can be processed similarly)
#
# See inst/exdata/mitomes/rCRS.fasta and inst/exdata/mitomes/rCRS.R 
#
if (FALSE) { 

  library(gmapR)
  library(rtracklayer)
  makeGmapGenome <- function(name, destDir="/home/tim/Dropbox/GmapGenomes") {
    fa <- system.file("extdata", paste0(name, ".fasta"),
                      package="ATACseeker", mustWork=TRUE) 
    fastaFile <- rtracklayer::FastaFile(fa)
    gmapGenome <- GmapGenome(fastaFile, create=TRUE)
    makeGmapGenomePackage(gmapGenome,
                          version="1.0.0", 
                          maintainer="<tim.triche@gmail.com>", 
                          author="Tim Triche, Jr.", 
                          destDir=destDir,
                          license="Artistic-2.0", 
                          pkgName=paste0("GmapGenome.Hsapiens.", name))
    return(gmapGenome)
  }

  # hg19 mitochondrial sequence: 
  GmapGenome.Hsapiens.hg19.chrM <- makeGmapGenome("hg19.chrM")

  # hg38 mitochondrial sequence: 
  GmapGenome.Hsapiens.hg38.chrM <- makeGmapGenome("hg38.chrM")

  # GRCh38 mitochondrial sequence:
  GmapGenome.Hsapiens.GRCh38.MT <- makeGmapGenome("GRCh38.MT")

}
