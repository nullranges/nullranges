#' Segmentation based on gene density
#'
#' @param x the input gene GRanges
#' @param n the number of states
#' @param Ls segment length
#' @param deny GRanges of deny region
#' @param type the type of segmentation, either "cbs" (which will
#' use DNAcopy to segment) or "hmm" (which will use RcppHMM).
#' The packages are not imported by nullranges, but must be installed
#' by the user
#' @param plot_origin plot the gene density of original gene Granges
#' @param boxplot boxplot of gene density for each state
#'
#' @import ggplot2
#' @importFrom plyranges filter join_overlap_intersect
#'
#' @export
segment_density <- function(x, n, Ls = 1e6, deny, type = c("cbs", "hmm"), plot_origin = TRUE, boxplot = FALSE) {
  query <- tileGenome(seqlengths(x)[seqnames(x)@values], tilewidth = Ls, cut.last.tile.in.chrom = TRUE)
  ## gap will create whole chromosome length ranges
  gap <- gaps(deny,end = seqlengths(x)) %>%
    plyranges::filter(strand=="*") 
  query_accept <- plyranges::join_overlap_intersect(query,gap2) %>% filter(width > Ls / 10)
  counts_nostand <- countOverlaps(query_accept, x)
  counts <- counts_nostand/width(query_accept2)*Ls
  eps <- rnorm(length(counts), 0, .2)

  if (plot_origin) {
    print(hist(counts, breaks = 50))
    # a<-seqnames(query)
    # b<-rep(a@values,a@lengths)
    print(plot(sqrt(counts) + eps))
  }

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
    scna <- DNAcopy::segment(cna,
      verbose = 1
    )
    seq <- with(scna$output, rep(seg.mean, num.mark))
    q <- quantile(seq, .95)
    seq2 <- pmin(seq, q)
    km <- kmeans(seq2, n)
    mcols(query_accept)$states <- km$cluster
    # TODO make optional to plot
    plot(sqrt(counts) + eps, col = km$cluster)
    
  } else if (type == "hmm") {

    if (!requireNamespace("RcppHMM", quietly=TRUE)) {
      stop("type='hmm' requires installing the Bioconductor package 'RcppHMM'")
    }
    
    hmm <- RcppHMM::initPHMM(n)
    hmm <- RcppHMM::learnEM(hmm,
      counts,
      iter = 400,
      delta = 1e-5,
      print = TRUE
    )
    v <- as.integer(factor(RcppHMM::viterbi(hmm, counts), levels = hmm$StateNames))
    plot(sqrt(counts) + eps, col = v)
    mcols(query_accept)$states <- v
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

  return(sort(seg))
}
