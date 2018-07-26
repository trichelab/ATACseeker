#' Generate bins for A/B compartment estimation
#'
#' Generate bins across a user defined chromosome for A/B compartment estimation.
#' A/B compartment estimation can be used for non-supported genomes if chr.end is set.
#' 
#' This function is used to generate a list object to be passed to getCorMatrix
#'
#' @param x      A p x n matrix where p (rows) = loci and n (columns) = samples/cells
#' @param genloc    GRanges object that contains corresponding genomic locations of the loci
#' @param chr      Chromosome to be analyzed
#' @param chr.start    Starting position (in bp) to be analyzed
#' @param chr.end    End position (in bp) to be analyzed
#' @param res    Binning resolution (in bp)
#' @param FUN    Function to be used to summarize information within a bin
#' @param genome    Genome corresponding to the input data ("hg19" or "mm10")
#' 
#' @return    A list object to pass to getCorMatrix
#' 
#' @import    GenomicRanges
#' @import    Homo.sapiens
#' @import    Mus.musculus
#' 
#' @export 

getBinMatrix <- function(x, genloc, chr = "chr1", chr.start = 0, chr.end = NULL, res = 100000, FUN=sum, genome = c("hg19", "mm10")){
  
  if (any(is.na(x))){
    stop("Matrix must not contain NAs")
  }
  if (nrow(x)!=length(genloc)){
    stop("Provided GRanges must have length equal to the matrix number of rows")
  }
  
  if (is.null(chr.end)) {
    if (genome == "hg19" | genome == "mm10") {
      genome <- match.arg(genome)
      chr.end <- switch(genome,
                        hg19 = seqlengths(Homo.sapiens)[chr],
                        mm10 = seqlengths(Mus.musculus)[chr])
    }
    else {
      message(paste0("Don't know what to do with ", genome))
      stop("If you'd like to use an unsupported genome, specify chr.end to an appropriate value...")
    }
  }
  start <- seq(chr.start, chr.end, res) #Build the possible bin ranges given resolution
  end <- c(start[-1], chr.end) - 1L #Set the end ranges for desired resolution
  
  #Build up the genomic ranges object given chr, start, end, and resolution
  gr.bin <- GRanges(seqnames = chr,
                    ranges = IRanges(start = start, end = end))
  
  #Identify overlaps between the user defined GRanges object (loci) and bins
  ids <- findOverlaps(genloc, gr.bin, select="first")
  
  #Get the number of bins overlapping loci
  n <- length(gr.bin)
  message(paste0(n, " bins created..."))
  
  #User defined function to summarize data in the bins
  #TODO: allow for bin matrices to be generated for all chrs
  x.bin <- apply(x, 2, function(x) {
    zvec <- rep(0, n) #Generate a vector of zeroes
    a <- tapply(x, INDEX=ids, FUN=FUN) #Summarize data
    zvec[as.numeric(names(a))] <- a
    zvec
  })
  
  colnames(x.bin) <- colnames(x) #Set colnames
  wh <- rowSums(x.bin) != 0 #Filter out empty bins
  
  #Subset the non-empty bins
  x.bin <- x.bin[wh,]
  gr.bin  <- gr.bin[wh]
  
  return(list(gr=gr.bin, x=x.bin))
}
