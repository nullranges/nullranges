#' Bootstrap GRanges
#'
#' @param x the input GRanges
#' @param type the type of null generation
#' @param L_b the length of the block
#' 
#' @export
bootstrap_granges <- function(x, type=c("bootstrap", "permute"), within_chrom=TRUE, L_b) {
  type <- match.arg(type)
  chrom_lens <- seqlengths(x)
  tab <- table(seqnames(x))
  chroms <- names(tab)
  if (within_chrom) {
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
    x_prime <- do.call(c, res)
  } else {
    # here we will send the ranges from multiple
    # chromosomes to a single long chromosome,
    # then perform block bootstrap / permutation,
    # then map back to the original chromosomes.
    
    # TODO: code assumes sorted 'x'
    stopifnot(all(x == GenomicRanges::sort(x)))

    # ranges on one long chromosome:
    r <- map_chroms_to_line(x)
    L_s <- sum(seqlengths(x))
    r_prime <- if (type == "bootstrap") {
                 bootstrap_iranges(r, L_s, L_b)
               } else {
                 permute_blocks_iranges(r, L_s, L_b)
               }
    r_prime <- GenomicRanges::sort(r_prime)
    # map the bootstrapped / permuted ranges back to chromosomes:
    x_prime <- map_line_to_chroms(r_prime, x)
  }
  x_prime
}

map_chroms_to_line <- function(x) {
  L_c <- unname(seqlengths(x))
  chrom_shift <- c(0,cumsum(L_c)[-length(L_c)])
  shift(ranges(x), rep(chrom_shift, seqnames(x)@lengths))
}

map_line_to_chroms <- function(r_prime, x) {
  L_c <- seqlengths(x)
  chrom_shift <- c(0,cumsum(unname(L_c))[-length(L_c)])
  chrom_blocks <- successiveIRanges(width=L_c)
  # need to deal with multiple overlaps...
  idx <- findOverlaps(r_prime, chrom_blocks, select="first")
  r_on_chroms <- shift(r_prime, -chrom_shift[idx])
  chroms <- seqlevels(x)[idx]
  keep <- start(r_on_chroms) >= 1 & end(r_on_chroms) <= L_c[chroms]
  GRanges(chroms[keep], r_on_chroms[keep], seqlengths=L_c)
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
