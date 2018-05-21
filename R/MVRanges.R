#' like a VRanges, but for mitochondria
#' 
#' @import VariantAnnotation
#' 
#' @exportClass MVRanges
setClass("MVRanges", 
         representation(coverage="numeric"),
         contains="VRanges")


#' wrap a VRanges for mitochondrial use
#'
#' @param   vr    the VRanges
#' @param   covg  estimated coverage
#'
#' @return        an MVRanges
#' 
#' @export
MVRanges <- function(vr, coverage) new("MVRanges", vr, coverage=coverage)


#' estimated read coverage (recorded from the called MAlignments)
#' 
#' @param x   an MVRanges
#' 
#' @return    estimated coverage (numeric) from the called MAlignments
#'
#' @export
setMethod("coverage", signature(x="MVRanges"), function(x) x@coverage)


#' types of variants (for feeding to MitImpact and friends) 
#' 
#' @param x   an MVRanges
#' 
#' @return    a vector of variant types ("SNV" or "indel")
#'
#' @export
setMethod("type", signature(x="MVRanges"), 
          function(x) ifelse(width(x) == 1, "SNV", "indel"))


#' positions of variants (for feeding to MitImpact and friends) 
#' 
#' @param x   an MVRanges
#' 
#' @return    a CHARACTER vector of variant position strings for MitImpact
#'
#' @export
setMethod("pos", signature(x="MVRanges"), 
          function(x) {
            gsub(paste0(seqlevels(x),":"), "", as.character(granges(x)))
          })


#' display variant calls with overall read coverage estimate
#'
#' @param x   an MVRanges
#' 
#' @export
setMethod("show", signature(object="MVRanges"),
          function(object) {
            callNextMethod()
            cat(paste0("  genome: ", unique(genome(object))))
            if ("anno" %in% names(metadata(object))) {
              cat(" (see metadata(object)$annotation)")
            }
            cat("\n")
            cat(paste0("  coverage: ~", round(coverage(object)), "x"), "\n")
          })

#' simple annotations of variants (using TxDb.Hsapiens.NCBI.rCRS)
#'
#' nb. Xiaowu's consequences are Synonymous, Missense, LoF, rRNA, tRNA, D-loop
#' so we kind of have to have D-loop annotation to tabulate these.  In order to
#' avoid doing this constantly, we cache the annotation in the MVRange's object
#' metadata() list. 
#' 
#' for reference, anno_rCRS was prepared as follows:
#'
#'   columns <- c("TXNAME","GENEID")
#'   anno <- transcripts(TxDb.Hsapiens.NCBI.rCRS, columns=columns)
#'   seqlevels(anno) <- "chrM" # just in case 
#'   anno$GENEID <- unlist(anno$GENEID) 
#'   is.na(anno[grep("id-", anno$TXNAME)]$GENEID) <- TRUE
#'   names(anno) <- ifelse(is.na(anno$GENEID), '', paste0("MT-", anno$GENEID))
#'
#'   # map gene names to expected symbols
#'   canon <- c("MT-COX1"="MT-CO1", "MT-COX2"="MT-CO2", "MT-COX3"="MT-CO3",
#'              "MT-CYTB"="MT-CYB")
#'   for (i in names(canon)) names(anno) <- sub(i, canon[i], names(anno))
#'
#'   # label the regions hit by various alterations 
#'   anno$region <- ifelse(!grepl("id-", anno$TXNAME), "coding",
#'                         ifelse(grepl("id-TRN", anno$TXNAME), "tRNA", "rRNA"))
#'   dLoop1 <- GRanges("chrM", IRanges(1, 576))
#'   seqinfo(dLoop1) <- seqinfo(anno)
#'   dLoop1$TXNAME <- NA_character_
#'   dLoop1$GENEID <- NA_character_
#'   dLoop1$region <- "D-loop"
#'   dLoop2 <- dLoop1
#'   end(dLoop2) <- 16569
#'   start(dLoop2) <- 16024
#'   anno <- c(dLoop1, anno, dLoop2)
#'
#' predictCoding(MVRanges, TxDb, BSgenome) still requires the TxDb and BSgenome 
#' objects for rCRS, so this just saves us a little time in regenerating them.
#' 
#' @param object  an MVRanges
#' 
#' @return        an MVRanges, with annotations from TxDb.Hsapiens.NCBI.rCRS
#'
#' @export
setMethod("annotation", signature(object="MVRanges"), 
          function(object) {

            if ("annotation" %in% names(metadata(object))) {
              anno <- metadata(object)$annotation 
            } else { 
              data(anno_rCRS)
              metadata(object)$annotation <- anno
            }
            
            if (unique(genome(object)) != "rCRS") {
              message("Lifting to rCRS...")
              object <- rCRS(object)
            }

            ol <- findOverlaps(object, anno)
            object$gene <- NA_character_
            object[queryHits(ol)]$gene <- names(anno)[subjectHits(ol)] 
            object$region <- NA_character_
            object[queryHits(ol)]$region <- anno[subjectHits(ol)]$region
            return(object)
          })
