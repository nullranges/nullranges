unseg_bootstrap <- function(x, L_b,
                            type = c("bootstrap", "permute"),
                            within_chrom = FALSE) {
  type <- match.arg(type)
  chrom_lens <- seqlengths(x)
  # TODO: what kind of solution do we want for this?
  # it's probably not a good idea to have L_b larger than one of the sequences
  stopifnot(all(chrom_lens >= L_b))
  tab <- table(seqnames(x))
  chroms <- names(tab)
  # first, simple case: sampling or permuting within chrom

  # TODO this loses the metadata columns right now...

  if (within_chrom) {
    r_primes <- lapply(chroms, function(chr) {
      # the length of the chromosomes
      L_s <- chrom_lens[[chr]]
      # the ranges on this chromosome
      r <- x[seqnames(x) == chr]
      r_prime <- if (type == "bootstrap") {
        unseg_bootstrap_within_chrom(r, L_b, L_s, chr)
      } else if (type == "permute") {
        unseg_permute_within_chrom(r, L_b, L_s, chr)
      }
      r_prime
    })
    x_prime <- do.call(c, r_primes)
  } else {
    # TODO: code assumes sorted 'x'
    stopifnot(all(x == GenomicRanges::sort(x)))
    # TODO: do we have to worry about missing seqlevels?
    # L_s is now a vector of the chromosome lengths
    L_s <- chrom_lens
    x_prime <- if (type == "bootstrap") {
      unseg_bootstrap_across_chrom(x, L_b, L_s)
    } else if (type == "permute") {
      unseg_permute_across_chrom(x, L_b, L_s)
    }
  }
  x_prime
}

# Block bootstrap GRanges within chromosome
#
# @param x the input GRanges
# @param L_b the length of the blocks
# @param L_s the length of the segment (chromosome)
# @param chr the name of the chromosome
unseg_bootstrap_within_chrom <- function(x, L_b, L_s=NULL, chr=NULL) {
  if (is.null(chr)) {
    chr <- as.character(seqnames(x)[1])
  }
  if (is.null(L_s)) {
    L_s <- seqlengths(x)[[chr]]
  }
  # blocks allowed to go over L_s
  n <- ceiling(L_s / L_b)
  # these blocks are the 'bait' to capture features in 'x'
  random_blocks <- IRanges::IRanges(start = round(runif(n, 1, (n - 1) * L_b + 1)),
                                    width = L_b)
  # where those blocks will move to
  rearranged_blocks <- IRanges::successiveIRanges(width(random_blocks))
  # the shift needed to move them
  block_shift <- start(rearranged_blocks) - start(random_blocks)
  # use the bait to sample features in 'x'
  fo <- IRanges::findOverlaps(random_blocks, ranges(x))
  # shift the ranges in those bait blocks
  suppressWarnings({
    x_prime <- shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
  })
  x_prime <- trim(x_prime)
  # sort outgoing ranges?
  sort(x_prime)
}

# Permute blocks of GRanges within chromosome
#
# @param x the input GRanges
# @param L_b the length of the blocks
# @param L_s the length of the segment (chromosome)
# @param chr the name of the chromosome
unseg_permute_within_chrom <- function(x, L_b, L_s=NULL, chr=NULL) {
  if (is.null(chr)) {
    chr <- as.character(seqnames(x)[1])
  }
  if (is.null(L_s)) {
    L_s <- seqlengths(x)[[chr]]
  }
  # blocks allowed to go over L_s
  n <- ceiling(L_s / L_b)
  # blocks tiling the chromosome
  blocks <- IRanges::successiveIRanges(rep(L_b, n))
  # which block do features in 'x' fall into?
  mcols(x)$block <- IRanges::findOverlaps(ranges(x), blocks, select = "first")
  # permutation order
  perm <- sample(n)
  # where those blocks will move to
  permuted_blocks <- IRanges::successiveIRanges(width(blocks))[perm]
  # the shift needed to move them
  block_shift <- start(permuted_blocks) - start(blocks)
  # shift the features in 'x'
  suppressWarnings({
    x_prime <- shift(x, block_shift[mcols(x)$block])
  })
  x_prime <- trim(x_prime)
  # sort outgoing ranges?
  GenomicRanges::sort(x_prime)
}

# Block bootstrap GRanges across chromosome
#
# @param x the input GRanges
# @param L_b the length of the blocks
# @param L_s the lengths of the chromosomes
unseg_bootstrap_across_chrom <- function(x, L_b, L_s) {
  # blocks allowed to go over L_s
  n_per_chrom  <- ceiling(L_s / L_b)
  # total number of tiling blocks
  n <- sum(n_per_chrom)
  # sample according to chromosome length
  random_chr <- sample(names(L_s), n, replace=TRUE, prob=L_s)
  # random starts within the tiled regions
  random_start <- round(runif(n, 1, (n_per_chrom[random_chr] - 1) * L_b + 1))
  # these blocks are the 'bait' for capturing features in 'x'
  random_blocks <- GenomicRanges::GRanges(seqnames=random_chr,
                           ranges=IRanges::IRanges(start = random_start, width = L_b))
  # where those blocks will move to
  rearranged_blocks <- GenomicRanges::tileGenome(seqlengths=L_s,
                                                 tilewidth=L_b,
                                                 cut.last.tile.in.chrom=TRUE)
  # use the bait to sample features in 'x'
  fo <- GenomicRanges::findOverlaps(random_blocks, x)
  # x has been sampled multiple times
  x_mult_hits <- x[subjectHits(fo)]
  # label which 'bait' block each feature hit in the re-sampling
  mcols(x_mult_hits)$block <- queryHits(fo)
  # shift the ranges in those bait blocks
  shift_and_swap_chrom(x_mult_hits, seqnames(rearranged_blocks),
                       start(random_blocks), start(rearranged_blocks))
}

# Permute blocks of GRanges across chomosome
#
# @param x the input GRanges
# @param L_b the length of the blocks
# @param L_s the lengths of the chromosomes
unseg_permute_across_chrom <- function(x, L_b, L_s) {
  blocks <- GenomicRanges::tileGenome(seqlengths=L_s, tilewidth=L_b,
                       cut.last.tile.in.chrom=TRUE)
    # pass along the full seqlengths of 'x'
  seqlengths(blocks) <- seqlengths(x)
  mcols(x)$block <- GenomicRanges::findOverlaps(x, blocks, select = "first")
  perm <- sample(length(blocks))
  rearranged_blocks <- blocks[perm]
  # this operation loses some ranges:
  # those that are mapped to permuted blocks that
  # are cut by `cut.last.tile.in.chrom`
  shift_and_swap_chrom(x, seqnames(rearranged_blocks),
                       start(blocks), start(rearranged_blocks))
}

# function moves featues in 'x' that fall into 'blocks'
# to the new block locations in 'rearranged_blocks', e.g.:
#
# 'x' in blocks --> 'x_prime' in rearranged_blocks
#
# and will also change the seqnames to 'chr_names'
shift_and_swap_chrom <- function(x, chr_names,
                                 random_blocks_start,
                                 rearranged_blocks_start) {
  stopifnot(length(rearranged_blocks_start) == length(random_blocks_start))
  block_shift <- rearranged_blocks_start - random_blocks_start
  idx <- mcols(x)$block
  chr_prime <- chr_names[idx]
  # this creates out-of-bound ranges
  # (but wait until we assign new chromosomes)
  suppressWarnings({
    x_prime <- shift(x, block_shift[idx])
    seqnames(x_prime) <- chr_prime
  })
  x_prime <- trim(x_prime)
  x_prime
}
