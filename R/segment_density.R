#' Genomic segmentation based on gene density
#'
#' @param x the input gene GRanges
#' @param n the number of states
#' @param L_s segment length
#' @param exclude GRanges of excluded region
#' @param type the type of segmentation, either "cbs" (which will
#' use DNAcopy to segment) or "hmm" (which will use RcppHMM).
#' The packages are not imported by nullranges, but must be installed
#' by the user
#' 
#' @return a GRanges with metadata columns containing:
#' \itemize{
#'   \item{state} {segmentation state}
#'   \item{counts} {average number of genes}
#' } 
#'
#' @importFrom plyranges filter join_overlap_intersect
#'
#' @export
segmentDensity <- function(x, n, L_s = 1e6, exclude,
                           type = c("cbs", "hmm")) {
  query <- GenomicRanges::tileGenome(seqlengths(x)[seqnames(x)@values],
    tilewidth = L_s,
    cut.last.tile.in.chrom = TRUE
  )
  # TODO: code assumes sorted 'x'
  if (any(x != GenomicRanges::sort(x))) {
    warning("unsorted x")
  }
  ## gap will create whole chromosome length ranges
  ## TODO: need to keep the gaps with same exclude strand,
  ## here is special case that all strand(exclude) ="*"
  gap <- gaps(exclude, end = seqlengths(x))
  gap <- plyranges::filter(gap, strand == "*")

  ## the region remove exclude regions
  query_accept <- filter(plyranges::join_overlap_intersect(query, gap),
                         width > L_s / 100)

  # TODO what does "nostand" mean? Answer: not standardized
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
  
  # recalculate counts on reduced segment region
  counts_nostand <- GenomicRanges::countOverlaps(seg, x, minoverlap = 8)
  counts <- counts_nostand / width(seg) * L_s
  mcols(seg)$counts <- counts
  return(seg)
}
