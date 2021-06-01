#' Genomic segmentation based on gene density
#'
#' @param x the input gene GRanges
#' @param seg the segmentation GRanges returned by \code{segmentDensity} function
#' @param L_s segment length
#' @param deny GRanges of deny region
#' @param type the type of plot returned. Either a segmentation plot 
#' included ranges information or a box plot. Default is both plots are displayed.
#' The y axis \code{"density"} represent square root of overlap counts within segment length. 
#' @param region GRanges of stricted region that want to be plotted. 
#'
#' @importFrom IRanges findOverlaps
#'
#' @export
plotSegment <- function(x, seg, L_s = 1e6, deny, type = c("ranges","boxplot"),
                        region = NULL) {
  
  counts_nostand <- GenomicRanges::countOverlaps(seg, x, minoverlap = 8)
  counts <- counts_nostand / width(seg) * L_s
  
  mcols(deny)$state <- "deny region"
  
  if(!is.null(region)){
    seg_fo <- findOverlaps(seg, region)
    seg <- seg[queryHits(seg_fo)]
    deny <- join_overlap_intersect(deny, region)
    counts <- counts[queryHits(seg_fo)]
  }
  full_query <- c(seg, deny)
  q <- quantile(sqrt(counts), .975)
  seq2 <- pmin(sqrt(counts), q)
  
  dat <- data.frame(
    chr = seqnames(full_query),
    counts = c(seq2, rep(0, length(deny))),
    state = factor(full_query$state),
    start = IRanges::start(full_query),
    end = IRanges::end(full_query))
  
  if("ranges" %in% type){
    p1 <- ggplot2::ggplot() +
      labs(col = "States", y = "density") +
      geom_segment(aes_string(x = "start", y = "counts", xend = "end", yend = "counts",col = "state"),
                   size = 5, data = dat, show.legend = TRUE) +
      theme_classic() 
    if(!is.null(region)){
      print(p1)
    }else{
      print(p1+facet_wrap(~chr, ncol = 4) +
              scale_x_continuous(
                name = "Location(Mb)",
                breaks = c(0.5e8, 1e8, 1.5e8, 2e8, 2.5e8),
                labels = c("50", "100", "150", "200", "250")
              ))
    }
    
  }
  if("boxplot" %in% type){
    states <- data.frame(state = factor(seg$state), count = seq2)
    p2 <- ggplot2::ggplot(aes_string(x = "state", y = "count"), data = states) +
      geom_boxplot() +
      theme_bw() +
      labs(x = "States", y = "density")
    print(p2)
  }
}