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
#' @return    a named CHARACTER vector of variant positions for MitImpact/HGVS
#'
#' @export
setMethod("pos", signature(x="MVRanges"), 
          function(x) {
            mtChr <- grep("(chrM|MT)", seqlevels(x), value=TRUE)
            loci <- gsub(paste0(mtChr,":"), "", as.character(granges(x)))
            return(loci)
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
            if ("annotation" %in% names(metadata(object))) {
              cat(" (try getAnnotations(object))")
            }
            cat(paste0(", ~", round(coverage(object)), "x read coverage")) 
            cat("\n")
          })


#' simple annotations of variants (using features from GenBank rCRS GFF record)
#'
#' nb. Xiaowu's consequences are Synonymous, Missense, LoF, rRNA, tRNA, D-loop
#' so we kind of have to have D-loop annotation to tabulate these.  In order to
#' avoid doing this constantly, we cache the annotation in the MVRange's object
#' metadata() list. This is used in the predictCoding() method for MVRanges.
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
#' @return        an MVRanges, with annotations from rCRS
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


#' simple helper to retrieve annotations (if present) from an MVRanges
#'
#' @param annotations   an MVRanges
#'
#' @return              a GRanges of annotations, if present, or else NULL
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


#' simple helper to retrieve PASS'ing, protein-coding variants from MVRanges
#'
#' @param x   an MVRanges
#'
#' @return    subset of MVRanges passing filters within protein-coding regions
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


#' predict coding mutations (annotating the most likely consequences)
#'
#' Given an rCRS-situated set of variant calls, predict which may be coding,
#' and which are synonymous (so as to compute dN/dS and annotate regionally).
#' 
#' Like the other predictCoding() methods, one can expect mcols() containing:
#' 
#' - varAllele   (reverse complemented if the subject is on the negative strand)
#' - QUERYID     (usually blank; rCRS liftOver step can make indices confusing)
#' - TXID        (usually blank; each mitochondrial gene has one transcript)
#' - CDSID       (usually blank; each mitochondrial gene has one known CDS) 
#' - GENEID      (usually blank; mitochondrial genes are just ID'ed as 1-13)
#' - CDSLOC      
#' - PROTEINLOC
#' - CONSEQUENCE (synonymous/nonsynonymous/frameshift/nonsense/not translated)
#' - REFCODON
#' - VARCODON
#' - REFAA
#' - VARAA
#'
#' Unlike the other predictCoding methods, many of these will usually be empty. 
#'
#' @param  query  an MVRanges
#'
#' @return        an rCRS MVRanges with CONSEQUENCE, MTCSQ, (REF/VAR)(CODON/AA)
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


#' summarize consensus pathogenicity estimates from MitImpact (2.0 or newer)
#'
#' Note: this method requires an internet connection. For now, at least.
#'
#' @param query   an MVRanges, almost always after predictCoding(MVRanges)
#'
#' @return        MitImpact output for variants where information is available 
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
