#' Combine nearby regions within same states for a user given input
#'
#' @param x the input gene GRanges
#' @param n the number of states and states are annoted as integer
#' @param col the name of the column pointing to segment states. Default is "state"
#'
#' @return a GRanges with metadata columns containing:
#' \itemize{
#'   \item{state} {segmentation state}
#' }
#'
#' @examples
#'
#' n <- 10000
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges(round(
#'   c(runif(n/4,1,991), runif(n/4,1001,3991),
#'     runif(n/4,4001,4991), runif(n/4,7001,9991))),
#'   width=10), seqlengths=c(chr1=10000))
#' gr$name = rep(1:4,each=10)
#' gr <- sort(gr)
#' seg <- reduceSegment(gr, n=4, col="name")
#'
#' @export
reduceSegment <- function(x, n, col="state") {

  seg <- do.call(c, lapply(seq_len(n), function(s) {
    x <- reduce(x[mcols(x)[,col] == s])
    mcols(x)$state <- s
    x
  }))

  seg <- sortSeqlevels(seg)
  seg <- GenomicRanges::sort(seg)

  return(seg)
}
