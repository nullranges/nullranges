#' Genomic segmentation based on gene density
#'
#' @param x the input gene GRanges
#' @param n the number of states
#' @param L_s segment length
#' @param deny GRanges of deny region
#' @param type the type of segmentation, either "cbs" (which will
#' use DNAcopy to segment) or "hmm" (which will use RcppHMM).
#' The packages are not imported by nullranges, but must be installed
#' by the user
#' @param returnPlot indicate whether return segmentation by chromosome
#' plot and boxplot. If true, plot element will be stored in the result list.
#'
#' @import ggplot2
#' @importFrom plyranges filter join_overlap_intersect
#'
#' @export
segmentDensity <- function(x, n, L_s = 1e6, deny,
                           type = c("cbs", "hmm"), returnPlot = FALSE) {
  query <- GenomicRanges::tileGenome(seqlengths(x)[seqnames(x)@values],
    tilewidth = L_s,
    cut.last.tile.in.chrom = TRUE
  )
  # TODO: code assumes sorted 'x'
  if (any(x != GenomicRanges::sort(x))) {
    warning("unsorted x")
  }
  ## gap will create whole chromosome length ranges
  ## TODO: need to keep the gaps with same deny strand,
  ## here is special case that all strand(deny) ="*"
  gap <- gaps(deny, end = seqlengths(x)) %>%
    plyranges::filter(strand == "*")

  ## the region remove deny regions
  query_accept <- plyranges::join_overlap_intersect(query, gap) %>%
    filter(width > L_s / 100)

  # TODO what does "nostand" mean?
  counts_nostand <- GenomicRanges::countOverlaps(query_accept, x, minoverlap = 8)
  counts <- counts_nostand / width(query_accept) * L_s

  if (type == "cbs") {
    if (!requireNamespace("DNAcopy", quietly = TRUE)) {
      stop("type='cbs' requires installing the Bioconductor package 'DNAcopy'")
    }

    cna <- DNAcopy::CNA(matrix(sqrt(counts), ncol = 1),
      chrom = as.character(seqnames(query_accept)),
      maploc = start(query_accept),
      data.type = "logratio",
      presorted = TRUE
    )
    smoothed.CNA.object <- DNAcopy::smooth.CNA(cna)
    scna <- DNAcopy::segment(smoothed.CNA.object,
      verbose = 1
    )

    seq <- rep(scna$output$seg.mean, scna$output$num.mark)
    km <- kmeans(seq, n)
    mcols(query_accept)$states <- km$cluster
  } else if (type == "hmm") {
    if (!requireNamespace("RcppHMM", quietly = TRUE)) {
      stop("type='hmm' requires installing the Bioconductor package 'RcppHMM'")
    }
    counts2 <- array(counts, dim = c(1, length(counts), 1))
    hmm <- RcppHMM::initGHMM(n)
    hmm <- RcppHMM::learnEM(hmm,
      counts2,
      iter = 400,
      delta = 1e-5,
      print = FALSE
    )
    v <- as.integer(factor(RcppHMM::viterbi(hmm, matrix(counts, nrow = 1)),
                           levels = hmm$StateNames))
    mcols(query_accept)$states <- v
  }

  # Combine nearby regions within same states
  seg <- do.call(c, lapply(1:n, function(s) {
    x <- reduce(query_accept[query_accept$states == s])
    mcols(x)$state <- s
    x
  }))

  seg <- sortSeqlevels(seg)
  seg <- GenomicRanges::sort(seg)

  if (returnPlot) {
    mcols(deny2)$states <- "deny region"
    full_query <- c(query_accept, deny2)
    q <- quantile(sqrt(counts), .975)
    seq2 <- pmin(sqrt(counts), q)
    dat <- data.frame(
      chr = seqnames(full_query), counts = c(seq2, rep(0, length(deny2))),
      state = factor(full_query$states, loc = IRanges::mid(full_query))
    )
    dat2 <- data.frame(chr = seqnames(deny2),
                       start = start(deny2),
                       end = end(deny2))
    p1 <- ggplot2::ggplot() +
      labs(col = "States", y = "sqrt(counts)") +
      geom_point(aes_string(x = "loc", y = "counts", col = "state"),
                 alpha = 0.5, data = dat) +
      geom_segment(aes_string(x = "start", y = 0, xend = "end", yend = 0),
                   size = 5, data = dat2, show.legend = FALSE) +
      facet_wrap(~chr, ncol = 4) +
      theme_bw() +
      scale_x_continuous(
        name = "Location(Mb)",
        breaks = c(0.5e8, 1e8, 1.5e8, 2e8, 2.5e8),
        labels = c("50", "100", "150", "200", "250")
      )

    states <- data.frame(state = factor(query_accept$states), count = seq2)
    p2 <- ggplot2::ggplot(aes_string(x = "state", y = "count"), data = states) +
      geom_boxplot() +
      theme_bw() +
      labs(x = "States", y = "sqrt(counts)")
    return(list(seg = seg, p1 = p1, p2 = p2))
  } else {
    return(seg)
  }
}
