#' Segmentation based on gene density
#'
#' @param x the input gene GRanges
#' @param n the number of states
#' @param L_s segment length
#' @param deny GRanges of deny region
#' @param type the type of segmentation, either "cbs" (which will
#' use DNAcopy to segment) or "hmm" (which will use RcppHMM).
#' The packages are not imported by nullranges, but must be installed
#' by the user
#' @param plot_segment plot the gene density segmentation by chromosome
#' @param boxplot boxplot of gene density for each state
#'
#' @import ggplot2
#' @importFrom plyranges filter join_overlap_intersect
#'
#' @export
segmentDensity <- function(x, n, L_s = 1e6, deny, type = c("cbs", "hmm"),
                            plot_segment = TRUE, boxplot = FALSE) {
  query <- tileGenome(seqlengths(x)[seqnames(x)@values],
                      tilewidth = L_s,
                      cut.last.tile.in.chrom = TRUE)
  ## gap will create whole chromosome length ranges
  ## To do: need to keep the gaps with same deny strand, here is special case that all strand(deny) ="*"
  gap <- gaps(deny,end = seqlengths(x)) %>%
    plyranges::filter(strand=="*")
  ## the region remove deny regions
  query_accept <- plyranges::join_overlap_intersect(query,gap) %>%
    filter(width > L_s / 100)
  counts_nostand <- countOverlaps(query_accept, x, minoverlap = 8)
  counts <- counts_nostand/width(query_accept) * L_s
  eps <- rnorm(length(counts), 0, .2)

  if (type == "cbs") {

    if (!requireNamespace("DNAcopy", quietly=TRUE)) {
      stop("type='cbs' requires installing the Bioconductor package 'DNAcopy'")
    }

    cna <- DNAcopy::CNA(matrix(sqrt(counts) + eps, ncol = 1),
      chrom = as.character(seqnames(query_accept)), 
      maploc = start(query_accept),
      data.type = "logratio",
      presorted = TRUE
    )
    smoothed.CNA.object <- DNAcopy::smooth.CNA(cna)
    scna <- DNAcopy::segment(smoothed.CNA.object,
      verbose = 1
    )
    seq <- with(scna$output, rep(seg.mean, num.mark))
    q <- quantile(seq, .95)
    seq2 <- pmin(seq, q)
    km <- kmeans(seq2, n)
    mcols(query_accept)$states <- km$cluster
  } else if (type == "hmm") {

    if (!requireNamespace("RcppHMM", quietly=TRUE)) {
      stop("type='hmm' requires installing the Bioconductor package 'RcppHMM'")
    }
    counts2 <- array(counts,dim = c(1,length(counts),1))
    hmm <- RcppHMM::initGHMM(n)
    hmm <- RcppHMM::learnEM(hmm,
      counts2,
      iter = 400,
      delta = 1e-5,
      print = FALSE
    )
    v <- as.integer(factor(RcppHMM::viterbi(hmm, matrix(counts,nrow=1)), levels = hmm$StateNames))
    mcols(query_accept)$states <- v
  }
  
  if (plot_segment) {
    mcols(deny)$states <- "black list"
    full_query <- c(query_accept,deny)
    q <- quantile(sqrt(counts) + eps, .95)
    seq2 <- pmin(sqrt(counts) + eps, q)
    dat <- data.frame(chr = seqnames(full_query), counts = c(seq2,rep(0,length(deny))), 
                      state = full_query$states, loc = c(end(full_query)))
    p <- ggplot2::ggplot(aes(x = loc, y = counts, col = factor(state)), data = dat) +
      geom_point() + labs(col = "States") +
      scale_x_continuous(name = "Location(Mb)",
                         breaks = c(0.5e8,1e8,1.5e8,2e8,2.5e8),
                         labels = c("50","100","150","200","250"))+
      facet_wrap(~chr,ncol = 4)+ theme_bw()+ coord_equal(ratio=6e6) 
    print(p)
  }
  
  if (boxplot) {
    states <- data.frame(state = query_accept$states, count = counts)
    p <- ggplot2::ggplot(aes(x = factor(state), y = counts), data = states) +
      geom_boxplot()
    print(p)
  }
  
  # Combine nearby regions within same states
  seg <- do.call(c, lapply(1:n, function(s) {
    x <- reduce(query_accept[query_accept$states == s])
    mcols(x)$state <- s
    x
  }))

  seg <- sortSeqlevels(seg)
  seg <- sort(seg)
}
