#' Bootstrap GRanges
#'
#' @param x the input GRanges
#' @param L_b the length of the block
#' @param type the type of null generation
#' @param within_chrom whether to re-sample (bootstrap) ranges
#' only wwithin each chromosome, or whether to allow
#' cross-chromosome re-sampling
#'
#' @export
bootstrap_granges <- function(x, L_b, type = c("bootstrap", "permute"), within_chrom = TRUE) {
  type <- match.arg(type)
  chrom_lens <- seqlengths(x)
  tab <- table(seqnames(x))
  chroms <- names(tab)
  # first, simple case: sampling or permuting within chrom
  if (within_chrom) {
    res <- lapply(chroms, function(chr) {
      # the length of the chromosomes
      L_s <- chrom_lens[[chr]]
      # the ranges on this chromosome
      r <- ranges(x[seqnames(x) == chr])
      r_prime <- if (type == "bootstrap") {
        bootstrap_iranges(r, L_s, L_b)
      } else {
        permute_blocks_iranges(r, L_s, L_b)
      }
      GRanges(seqnames = chr, ranges = r_prime, seqlengths = chrom_lens)
    })
    x_prime <- do.call(c, res)
  } else {
    # TODO: code assumes sorted 'x'
    stopifnot(all(x == GenomicRanges::sort(x)))

    # TODO: code right now assumes these are the same... fix later
    stopifnot(all(sort(unique(seqnames(x))) == sort(seqlevels(x))))

    # L_s is now a vector of the chromosome lengths
    L_s <- chrom_lens
    x_prime <- if (type == "bootstrap") {
      block_bootstrap_granges(x, L_s, L_b)
    } else {
      permute_blocks_granges(x, L_s, L_b)
    }
    # sort outgoing ranges?
    x_prime <- GenomicRanges::sort(x_prime)
  }
  x_prime
}

#' Block bootstrap IRanges
#'
#' @param x the input ranges
#' @param L_s the length of the segment
#' @param L_b the length of the blocks
#'n
#' @export
bootstrap_iranges <- function(x, L_s, L_b) {
  # blocks allowed to go over L_s
  n <- ceiling(L_s / L_b)
  # these blocks are the 'bait' to capture features in 'x'
  random_blocks <- IRanges(start = round(runif(n, 1, (n - 1) * L_b + 1)), width = L_b)
  # where those blocks will move to
  rearranged_blocks <- successiveIRanges(width(random_blocks))
  # the shift needed to move them
  block_shift <- start(rearranged_blocks) - start(random_blocks)
  # use the bait to sample features in 'x'
  fo <- findOverlaps(random_blocks, x)
  # shift the ranges in those bait blocks
  x_prime <- shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
  # here trim any features that are not within [1, L_s]
  x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_s]
  # sort outgoing ranges?
  sort(x_prime)
}

#' Permute blocks IRanges
#'
#' @param x the input ranges
#' @param L_s the length of the segment
#' @param L_b the length of the blocks
#'
#' @export
permute_blocks_iranges <- function(x, L_s, L_b) {
  # blocks allowed to go over L_s
  n <- ceiling(L_s / L_b)
  # blocks tiling the chromosome
  blocks <- successiveIRanges(rep(L_b, n))
  # which block do features in 'x' fall into?
  mcols(x)$block <- findOverlaps(x, blocks, select = "first")
  # permutation order
  perm <- sample(n)
  # where those blocks will move to
  permuted_blocks <- successiveIRanges(width(blocks))[perm]
  # the shift needed to move them
  block_shift <- start(permuted_blocks) - start(blocks)
  # shift the features in 'x'
  x_prime <- shift(x, block_shift[mcols(x)$block])
  # here trim any features that are not within [1, L_s]
  x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_s]
  # sort outgoing ranges?
  sort(x_prime)
}

# Block bootstrap GRanges
#
# @param x the input ranges
# @param L_s the lengths of the chromosomes
# @param L_b the length of the blocks
#
block_bootstrap_granges <- function(x, L_s, L_b) {
  # blocks allowed to go over L_s
  n_per_chrom  <- ceiling(L_s / L_b)
  # total number of tiling blocks
  n <- sum(n_per_chrom)
  # sample according to chromosome length
  random_chr <- sample(names(L_s), n, replace=TRUE, prob=L_s)
  # random starts within the tiled regions
  random_start <- round(runif(n, 1, (n_per_chrom[random_chr] - 1) * L_b + 1))
  # these blocks are the 'bait' for capturing features in 'x'
  random_blocks <- GRanges(seqnames=random_chr,
                           ranges=IRanges(start = random_start, width = L_b))
  # where those blocks will move to
  rearranged_blocks <- tileGenome(seqlengths=L_s, tilewidth=L_b, cut.last.tile.in.chrom=TRUE)
  # use the bait to sample features in 'x'
  fo <- findOverlaps(random_blocks, x)
  # x has been sampled multiple times
  x_mult_hits <- x[subjectHits(fo)]
  # label which 'bait' block each feature hit in the re-sampling
  mcols(x_mult_hits)$block <- queryHits(fo)
  # shift the ranges in those bait blocks
  shift_and_swap_chrom(x_mult_hits, random_blocks, rearranged_blocks)
}

# Permute blocks GRanges
#
# @param x the input ranges
# @param L_s the lengths of the chromosomes
# @param L_b the length of the blocks
# 
permute_blocks_granges <- function(x, L_s, L_b) {
  blocks <- tileGenome(seqlengths=L_s, tilewidth=L_b, cut.last.tile.in.chrom=TRUE)
  mcols(x)$block <- findOverlaps(x, blocks, select = "first")
  perm <- sample(length(blocks))
  rearranged_blocks <- blocks[perm]
  # this operation loses some ranges:
  # those that are mapped to permuted blocks that
  # are cut by `cut.last.tile.in.chrom`
  shift_and_swap_chrom(x, blocks, rearranged_blocks)
}

# function moves featues in 'x' that fall into 'blocks'
# to the new block locations in 'rearranged_blocks', e.g.:
#
# 'x' in blocks --> 'x_prime' in rearranged_blocks
#
# and will also change the seqnames
shift_and_swap_chrom <- function(x, blocks, rearranged_blocks) {
  block_shift <- start(rearranged_blocks) - start(blocks)
  idx <- mcols(x)$block
  chr_prime <- seqnames(rearranged_blocks)[idx]
  # this creates out-of-bound ranges
  # (but wait until we assign new chromosomes)
  suppressWarnings({
    x_prime <- shift(x, block_shift[idx])
    seqnames(x_prime) <- chr_prime
  })
  x_prime <- trim(x_prime)
  x_prime <- x_prime[width(x_prime) > 0]
  x_prime
}
