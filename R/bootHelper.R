#' Segmentation based on one region
#'
#' This function makes a segmentation (GRanges) based on one
#' region of one chromosome (seqnames).
#'
#' @param x a single region as GRanges object
#' @param seqlength optional, the length of the chromosome,
#' if not provided, the function will attempt to pull this
#' using \code{genome(x)} and the \code{Seqinfo} function
#'
#' @return a segmentation (GRanges object) with the region of
#' interest designated as state 2, and the rest of the chromosome
#' as state 1.
#'
#' @examples
#'
#' library(GenomicRanges)
#' library(GenomeInfoDb)
#' x <- GRanges("chr1", IRanges(10e6+1,width=1e6))
#' genome(x) <- "hg19"
#' seg <- oneRegionSegment(x)
#' 
#' @export
oneRegionSegment <- function(x, seqlength) {
  stopifnot(length(x) == 1)
  chrom <- as.character(seqnames(x))
  g <- genome(x)[[chrom]]
  if (missing(seqlength)) {
    if (!is.na(g)) {
      si <- GenomeInfoDb::Seqinfo(genome=g)
      seqlength <- seqlengths(si[chrom])[[1]]
    } else {
      stop("function requires seqlength argument or Seqinfo-supported genome")
    }
  }
  GRanges(seqnames(x),
          IRanges(c(1,start(x),end(x)+1),
                  c(start(x)-1,end(x),seqlength)),
          state=c(1,2,1))
}
