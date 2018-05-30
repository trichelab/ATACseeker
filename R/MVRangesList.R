#' like a VRangesList, but for mitochondria
#' 
#' @import VariantAnnotation
#' @import Biostrings
#' @import S4Vectors
#' @import chromVAR
#' 
#' @exportClass MVRangesList
setClass("MVRangesList", contains="SimpleVRangesList")


#' wrap a VRangesList for mitochondrial use
#'
#' @param ...     the MVRanges elements forming the MVRangesList
#'
#' @return        the MVRangesList
#' 
#' @export
MVRangesList <- function(...) {
  new("MVRangesList", GenomicRangesList(...), elementType = "MVRanges")
}


#' MVRangesList methods (centralized).
#'
#' `counts`               returns fragment counts, if any
#' `counts<-`             adds or updates fragment counts
#' `coverage`             returns estimated coverage for each element
#' `encoding`             returns mutations in coding regions for each element
#' `granges`              returns mildly annotated aggregates of variant sites
#' `summarizeVariants`    attempts mass functional annotation of variant sites
#' 
#' @param x             an MVRangesList (for some methods)
#' @param object        an MVRangesList (for other methods)
#' @param value         a RangedSummarizedExperiment with matching colnames
#' @param annotations   a RangedSummarizedExperiment with motif count matches
#' @param filterLowQual optional argument to `granges` and `summarizeVariants`
#'
#' @name  MVRangesList-methods
NULL


#' @rdname    MVRangesList-methods
#' @export
setMethod("coverage", signature(x="MVRangesList"), 
          function(x) sapply(x, coverage))


#' @rdname    MVRangesList-methods
#' @export
setReplaceMethod("counts", 
                 signature(object="MVRangesList", 
                           value="RangedSummarizedExperiment"),
                 function(object, value) {
                   if (!identical(names(object), colnames(value))) {
                     stop("Error: colnames(value) doesn't match names(object)!")
                   } else if (!"counts" %in% names(assays(value))) {
                     stop("Error: value must have an assay named `counts`!")
                   } else {
                     columns <- names(object)
                     metadata(object)$counts <- filterPeaks(value[, columns])
                     return(object)
                   }
                 })


#' @rdname    MVRangesList-methods
#' @export
setMethod("counts", signature(object="MVRangesList"), 
          # it turns out that filtering may be needed on egress:
          function(object) filterPeaks(metadata(object)$counts))


#' @rdname    MVRangesList-methods
#' @export
setMethod("encoding", signature(x="MVRangesList"), 
          function(x) MVRangesList(lapply(x, encoding)))


#' @rdname    MVRangesList-methods
#' @export
setMethod("show", signature(object="MVRangesList"),
          function(object) {
            callNextMethod()
            coverages <- paste0(round(unname(sapply(object, coverage))), "x")
            cat(S4Vectors:::labeledLine("coverage", coverages))
            if ("counts" %in% names(metadata(object))) {
              peaks <- nrow(metadata(object)$counts)
              cat(ifelse("bias" %in% names(rowData(counts(object))),
                  "Bias-corrected ", "Raw "))
              cat("fragment counts at", peaks, "peaks are available from",
                  "counts(object).\n")
            }
          })


#' @rdname    MVRangesList-methods
#' @export
setMethod("granges", signature(x="MVRangesList"),
          function(x, filterLowQual=TRUE) {
            data(hg19TorCRS)
            data(chrominfo.rCRS)
            # pull in annotations
            anno <- suppressMessages(getAnnotations(annotation(x[[1]]))) 
            if (filterLowQual == TRUE) {
              message("Filtering out low-quality calls...")
              x <- MVRangesList(sapply(x, subset, PASS))
            }
            message("Aggregating variants...")
            gr <- unlist(GRangesList(sapply(x, granges)))
            seqlevelsStyle(gr) <- "UCSC" # chrM
            mtGenome <- unique(genome(gr))
            if (mtGenome %in% c("rCRS","GRCh38","hg38")) {
              seqinfo(gr) <- chrominfo.rCRS # identical save for name 
            } else if (mtGenome == "hg19") {
              message("Lifting variants to rCRS...")
              gr <- sort(unlist(liftOver(gr, hg19TorCRS)))
              seqinfo(gr) <- chrominfo.rCRS # as with TxDB
            } else if (mtGenome != "rCRS") { 
              stop("Unsupported genome: ", mtGenome)
            }
            gr <- reduce(gr)
            annoGenome <- unique(genome(anno))
            newMtGenome <- unique(genome(gr))
            stopifnot(newMtGenome == annoGenome)
            metadata(gr)$annotation <- anno
            ol <- findOverlaps(gr, anno)
            message("Annotating variants by region...")
            gr$gene <- NA_character_
            gr[queryHits(ol)]$gene <- names(anno)[subjectHits(ol)] 
            gr$region <- NA_character_
            gr[queryHits(ol)]$region <- anno[subjectHits(ol)]$region
            return(gr)
          })


#' @rdname    MVRangesList-methods
#' @export
setMethod("summarizeVariants", 
          signature(query="MVRangesList","missing","missing"),
          function(query, ...) {
            
            # code duplication! refactor
            getRangedImpact <- function(pos) {
              url <- paste("http://mitimpact.css-mendel.it", "api", "v2.0",
                           "genomic_position", pos, sep="/")
              res <- as.data.frame(read_json(url, simplifyVector=TRUE)$variants)
              if (nrow(res) > 0) {
                res$genomic <- with(res, paste0("g.", Start, Ref, ">", Alt))
                res$protein <- with(res, paste0("p.",AA_ref,AA_position,AA_alt))
                res$Consequence <- with(res, 
                                        paste(Gene_symbol, 
                                              paste0(AA_ref,
                                                     AA_position,
                                                     AA_alt)))
                res[, c("genomic","protein","Start",
                        "Ref","Alt","Codon_substitution","dbSNP_150_id",
                        "Mitomap_Phenotype","Mitomap_Status",
                        "Gene_symbol","OXPHOS_complex",
                        "Consequence","APOGEE_boost_consensus","MtoolBox")]
              } else {
                return(NULL)
              }
            }

            gr <- granges(query, ...)
            names(gr) <- as.character(gr)
            message("Retrieving functional annotations for variants...")
            hits <- lapply(as.character(ranges(gr)), getRangedImpact)
            rsv <- do.call(rbind, hits[which(sapply(hits, length) > 0)])
            names(rsv) <- sub("Start", "start", names(rsv)) # grrr
            rsv$chrom <- "chrM"
            rsv$end <- rsv$start # FIXME
            res <- makeGRangesFromDataFrame(rsv, keep=TRUE)
            seqinfo(res) <- seqinfo(gr)
            return(res)

          })

