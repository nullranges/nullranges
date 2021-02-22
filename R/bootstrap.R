#' Bootstrap GRanges
#'
#' @param x the input GRanges
#' @param type the type of null generation
#' @param L_b the length of the block
#' 
#' @export
bootstrap_granges <- function(x, type=c("bootstrap", "permute"), L_b) {
  type <- match.arg(type)
  chrom_lens <- seqlengths(x)
  tab <- table(seqnames(x))
  chroms <- names(tab)
  res <- lapply(chroms, function(chr) {
    L_s <- chrom_lens[[chr]]
    r <- ranges(x[seqnames(x) == chr])
    r_prime <- if (type == "bootstrap") {
      bootstrap_iranges(r, L_s, L_b)
    } else {
      permute_blocks_iranges(r, L_s, L_b)
    }
    GRanges(seqnames=chr, ranges=r_prime, seqlengths=chrom_lens)
  })
  do.call(c, res)
}

#' Block bootstrap IRanges
#'
#' @param x the input ranges
#' @param L_s the length of the segment
#' @param L_b the length of the blocks
#'
#' @export
bootstrap_iranges <- function(x, L_s, L_b) {
  # first attempt: blocks go over L_s, and we will trim them later
  n <- ceiling(L_s/L_b) 
  random_blocks <- IRanges(start=round(runif(n, 1, (n-1)*L_b+1)), width=L_b)
  rearranged_blocks <- successiveIRanges(width(random_blocks))
  block_shift <- start(rearranged_blocks) - start(random_blocks)
  fo <- findOverlaps(random_blocks, x)
  x_prime <- shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
  x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_s]
  x_prime
}

#' Permute blocks IRanges
#'
#' @param x the input ranges
#' @param L_s the length of the segment
#' @param L_b the length of the blocks
#'
#' @export
permute_blocks_iranges <- function(x, L_s, L_b) {
  n <- ceiling(L_s/L_b)
  blocks <- successiveIRanges(rep(L_b, n))
  mcols(x)$block <- findOverlaps(x, blocks, select="first")
  perm <- sample(n)
  permuted_blocks <- successiveIRanges(width(blocks)[perm])
  permuted_blocks[perm] <- permuted_blocks
  block_shift <- start(permuted_blocks) - start(blocks)
  x_prime <- shift(x, block_shift[mcols(x)$block])
  x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_s]
  x_prime
}
