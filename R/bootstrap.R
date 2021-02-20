#' Bootstrap GRanges
#'
#' @export
bootstrap_granges <- function(x, type=c("bootstrap", "permute"), W) {
  type <- match.arg(type)
  Ls <- seqlengths(x)
  tab <- table(seqnames(x))
  chroms <- names(tab)
  res <- lapply(chroms, function(chr) {
    L <- Ls[[chr]]
    r <- ranges(x[seqnames(x) == chr])
    r_prime <- if (type == "bootstrap") {
      bootstrap_iranges(r, L, W)
    } else {
      permute_blocks_iranges(r, L, W)
    }
    GRanges(seqnames=chr, ranges=r_prime, seqlengths=Ls)
  })
  do.call(c, res)
}

#' Block bootstrap IRanges
#'
#' @export
bootstrap_iranges <- function(x, L, W) {
  # first attempt: blocks go over L, and we will trim them later
  n <- ceiling(L/W) 
  random_blocks <- IRanges(start=round(runif(n, 1, (n-1)*W+1)), width=W)
  rearranged_blocks <- successiveIRanges(width(random_blocks))
  block_shift <- start(rearranged_blocks) - start(random_blocks)
  fo <- findOverlaps(random_blocks, x)
  x_prime <- shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
  x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L]
  x_prime
}

#' Permute blocks IRanges
#'
#' @export
permute_blocks_iranges <- function(x, L, W) {
  n <- ceiling(L/W)
  blocks <- successiveIRanges(rep(W, n))
  mcols(x)$block <- findOverlaps(x, blocks, select="first")
  perm <- sample(n)
  permuted_blocks <- successiveIRanges(width(blocks)[perm])
  permuted_blocks[perm] <- permuted_blocks
  block_shift <- start(permuted_blocks) - start(blocks)
  x_prime <- shift(x, block_shift[mcols(x)$block])
  x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L]
  x_prime
}
