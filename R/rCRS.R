#' liftOver and simplify an hg19 mitochondrial genome (variants) to rCRS 
#'
#' @param mr    an MRanges of variant calls
#'
#' @return      an MRanges of variant calls against rCRS 
#' 
#' @import rtracklayer
#'
#' @export
rCRS <- function(mr) { 
  mtGenome <- unique(genome(mr))
  if (mtGenome %in% c("GRCh38","hg38")) {
    message("This MRanges is already against rCRS.")
    return(mr)
  } else { 
    data(hg19TorCRS) 
    mr <- sort(unlist(liftOver(mr, hg19TorCRS)))
    names(mr) <- paste0(as.character(mr), ":", ref(mr), ">", alt(mr))
    return(mr)
  }
}
