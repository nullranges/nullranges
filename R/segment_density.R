#' Segmentation based on gene density
#'
#' @param x the input gene GRanges
#' @param n the number of states
#' @param Ls segment length
#' @param deny GRanges of deny region
#' @param type the type of segmentation
#' @param plot_origin plot the gene density of original gene Granges
#' @param boxplot boxplot of gene density for each state
#'
#' @import DNAcopy
#' @import RcppHMM
#' @import ggplot2
#' @importFrom plyranges filter join_overlap_intersect
#'
#' @export
segment_density <- function(x, n, Ls = 1e6, deny, type = c("CBS", "HMM"), plot_origin = TRUE, boxplot = FALSE) {
  query <- tileGenome(seqlengths(x)[seqnames(x)@values], tilewidth = Ls, cut.last.tile.in.chrom = TRUE)
  gap <- gaps(deny,end = seqlengths(x)) %>% plyranges::filter(strand=="*") ## gap will create whole chromosome length ranges
  query_accept <- plyranges::join_overlap_intersect(query,gap)
  counts_nostand <- countOverlaps(query_accept, x)
  counts <- counts_nostand/width(query_accept)*Ls
  # counts <- log2(counts+1)
  eps <- rnorm(length(counts), 0, .2)
  if (plot_origin) {
    print(hist(counts, breaks = 50))
    # a<-seqnames(query)
    # b<-rep(a@values,a@lengths)
    print(plot(sqrt(counts) + eps))
  }

  if (type == "CBS") {
    cna <- CNA(matrix(sqrt(counts) + eps, ncol = 1),
      chrom = as.character(seqnames(query_accept)), # wont work for X,Y,MT
      maploc = start(query_accept),
      data.type = "logratio",
      presorted = TRUE
    )
    scna <- segment(cna,
      undo.splits = "sdundo",
      undo.SD = 1.5,
      verbose = 1
    )
    seq <- with(scna$output, rep(seg.mean, num.mark))
    # plot(scna)
    # plot(seq)
    q <- quantile(seq, .95)
    seq2 <- pmin(seq, q)
    # plot(seq2)
    km <- kmeans(seq2, n)
    query_accept$states <- km$cluster
    plot(sqrt(counts) + eps, col = km$cluster)
  } else {
    hmm <- initPHMM(n)
    hmm <- learnEM(hmm,
      counts,
      iter = 400,
      delta = 1e-5,
      print = TRUE
    )
    hmm
    v <- as.integer(factor(viterbi(hmm, counts), levels = hmm$StateNames))
    plot(sqrt(counts) + eps, col = v)
    query_accept$states <- v
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
