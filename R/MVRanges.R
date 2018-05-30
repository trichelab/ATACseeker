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

#' MVRanges methods (centralized).
#'
#' @name      MVRanges-methods
NULL


#' @rdname    MVRanges-methods
#' 
#' @param x   an MVRanges
#' 
#' @return    estimated coverage (numeric) from the called MAlignments
#'
#' @export
setMethod("coverage", signature(x="MVRanges"), function(x) x@coverage)


#' @rdname    MVRanges-methods
#' 
#' @param x   an MVRanges
#' 
#' @return    a vector of variant types ("SNV" or "indel")
#'
#' @export
setMethod("type", signature(x="MVRanges"), 
          function(x) ifelse(width(x) == 1, "SNV", "indel"))


#' @rdname    MVRanges-methods
#' 
#' @param x   an MVRanges
#' 
#' @return    a named CHARACTER vector of variant positions for MitImpact/HGVS
#'
#' @export
setMethod("pos", signature(x="MVRanges"), 
          function(x) {
            mtChr <- grep("(chrM|MT)", seqlevels(x), value=TRUE)
            loci <- gsub(paste0(mtChr,":"), "", as.character(granges(x)))
            return(loci)
          })


#' @rdname    MVRanges-methods
#'
#' @param x   an MVRanges
#' 
#' @export
setMethod("show", signature(object="MVRanges"),
          function(object) {
            callNextMethod()
            cat(paste0("  genome: ", unique(genome(object))))
            if ("annotation" %in% names(metadata(object))) {
              cat(" (try getAnnotations(object))")
            }
            cat(paste0(", ~", round(coverage(object)), "x read coverage")) 
            cat("\n")
          })


#' @rdname    MVRanges-methods
#' 
#' @param object  an MVRanges
#'
#' @export
setMethod("annotation", signature(object="MVRanges"), 
          function(object) {

            if (!"annotation" %in% names(metadata(object))) {
              if (unique(genome(object)) != "rCRS") {
                message("Lifting to rCRS...")
                object <- rCRS(object)
              }
              data(anno_rCRS)
              metadata(object)$annotation <- anno_rCRS
            }

            anno <- getAnnotations(object)
            ol <- findOverlaps(object, anno)
            object$gene <- NA_character_
            object[queryHits(ol)]$gene <- names(anno)[subjectHits(ol)] 
            object$region <- NA_character_
            object[queryHits(ol)]$region <- anno[subjectHits(ol)]$region
            return(object)

          })


#' @rdname    MVRanges-methods
#'
#' @param annotations   an MVRanges
#'
#' @export
setMethod("getAnnotations", signature(annotations="MVRanges"), 
          function(annotations) {
            anno <- metadata(annotations)$annotation
            if (is.null(anno)) {
              message("Unannotated! Try getAnnotations(annotation(object))).")
            }
            return(anno)
          })


#' @rdname    MVRanges-methods
#'
#' @param x   an MVRanges
#'
#' @import    Biostrings
#'
#' @export
setMethod("encoding", signature(x="MVRanges"), 
          function(x) {

            # limit the search 
            x <- annotation(x) # ensure it's rCRS
            x <- subset(x, PASS & region == "coding") 

            # fix issues
            data(rCRSeq)
            comp <- data.frame(ref=as.character(getSeq(rCRSeq, x)), alt=alt(x))
            keep <- with(comp, which(ref != alt))
            
            # return subset
            return(x[keep])

          })


#' @rdname    MVRanges-methods
#'
#' @param  query  an MVRanges
#'
#' @import        Biostrings
#'
#' @export
setMethod("predictCoding", # mitochondrial annotations kept internally
          signature(query="MVRanges", "missing", "missing", "missing"), 
          function(query, ...) {

            # setup:
            data(rCRSeq)
            query <- encoding(rCRS(query)) 
            MT_CODE <- getGeneticCode("SGC1")
            mtGenes <- subset(metadata(query)$annotation, region == "coding")
            mcols(mtGenes) <- DataFrame(DNA=getSeq(rCRSeq, mtGenes))
            mtGenes$AA <- translate(mtGenes$DNA, MT_CODE)
            ol <- findOverlaps(query, mtGenes)

            # execution:
            result <- granges(query)
            result$varAllele <- alt(query)

            stop("predictCoding(MVRanges) is not finished") 

          })


#' @rdname    MVRanges-methods
#'
#' @param query   an MVRanges, almost always after predictCoding(MVRanges)
#'
#' @import        jsonlite
#'
#' @export
setMethod("summarizeVariants", signature(query="MVRanges","missing","missing"),
          function(query, ...) {
           
            getImpact <- function(pos) {
              url <- paste("http://mitimpact.css-mendel.it", "api", "v2.0",
                           "genomic_position", sub("^g\\.", "", pos), sep="/")
              res <- as.data.frame(read_json(url, simplifyVector=TRUE)$variants)
              if (nrow(res) > 0) {
                res$genomic <- with(res, paste0("g.", Start, Ref, ">", Alt))
                res$protein <- with(res, paste0("p.",AA_ref,AA_position,AA_alt))
                res$change <- with(res, paste(Gene_symbol, protein))
                res[, c("genomic","protein","APOGEE_boost_consensus","MtoolBox",
                        "Mitomap_Phenotype","Mitomap_Status","OXPHOS_complex",
                        "dbSNP_150_id","Codon_substitution")]
              } else {
                return(NULL)
              }
            }

            hits <- lapply(pos(encoding(query)), getImpact)
            hits <- hits[which(sapply(hits, length) > 0)] 

            # be precise, if possible
            for (h in names(hits)) {
              j <- hits[[h]]
              if (h %in% j$genomic) hits[[h]] <- j[which(j$genomic == h),]
            }

            do.call(rbind, hits)

          })
