#' liftOver and simplify an hg19 mitochondrial genome (variants) to rCRS 
#'
#' @param mvr   an MVRanges of variant calls
#'
#' @return      an MVRanges of variant calls against rCRS 
#' 
#' @import rtracklayer
#'
#' @export
rCRS <- function(mvr) { 
  mtGenome <- unique(genome(mvr))
  if (mtGenome %in% c("GRCh38","hg38")) {
    message("This MVRanges is already against rCRS.")
    return(mvr)
  } else { 
    data(hg19TorCRS) 
    mvr <- sort(unlist(liftOver(mvr, hg19TorCRS)))
    names(mvr) <- paste0(as.character(mvr), ":", ref(mvr), ">", alt(mvr))
    # some variants called against hg19 chrM are really just the rCRS itself 
    message("Warning: 'cleanup' of hg19-to-rCRS artifactual 'SNPs' is not done")
    return(mvr)
  }
}
