#' Segmentated bootstrap GRanges
#'
#' @param seg segmentation GRanges a column("state") indicate segmentation states
#' @param x the input GRanges
#' @param L_c the length of the chromosome
#' @param L_s a vector of the length of each segmentation region
#' @param L_b the length of the block
#' 
#' @export
seg_bootstrap_iranges <- function(seg,x,L_c, L_s,L_b) {
  ns<-length(L_s)
  L_b0<-round(L_s*L_b/L_c) #block width for each segmentation state
  n <- ceiling(L_c/L_b) 
  r<-ranges(x)
  obj<-lapply(1:ns,function(m){ #loop over segmentation state
    L_bs<-L_b0[m] 
    seg2<-ranges(seg[seg$state == m])
    
    p<-seg2@width/sum(seg2@width)
    times<-round(n*p) # number of block within every piece segmentation
    
    start<-mapply(function(time,x, y) runif(time, x, y), times, seg2@start,end(seg2))
    if (is.list(start)) {
      start <- do.call(c, start)
    }
    start<-sample(start)
    random_blocks <- IRanges(start=start,width=L_bs)
    
    start_order<-mapply(function(x,y) seq(x,y,L_bs),seg2@start,end(seg2))
    if (is.list(start_order)) {
      start_order <- do.call(c, start_order)
    }
    rearranged_blocks <- IRanges(start=start_order,width=L_bs)
    
    block_shift <- start(rearranged_blocks) - start(random_blocks)
    
    fo <- findOverlaps(random_blocks, r)
    x_prime <- shift(r[subjectHits(fo)], block_shift[queryHits(fo)])
    x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_c]
    x_prime
  })
  return(obj)
}
