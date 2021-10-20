unseg_bootstrap <- function(y, L_b,
                            type = c("bootstrap", "permute"),
                            within_chrom = FALSE) {
  type <- match.arg(type)
  chrom_lens <- seqlengths(y)
  # TODO: what kind of solution do we want for this?
  # it's probably not a good idea to have L_b larger than one of the sequences
  stopifnot(all(chrom_lens >= L_b))
  tab <- table(seqnames(y))
  chroms <- names(tab)
  # first, simple case: sampling or permuting within chrom

  # TODO this loses the metadata columns right now...

  if (within_chrom) {
    r_primes <- lapply(chroms, function(chr) {
      # the length of the chromosomes
      L_s <- chrom_lens[[chr]]
      # the ranges on this chromosome
      r <- y[seqnames(y) == chr]
      r_prime <- if (type == "bootstrap") {
        unseg_bootstrap_within_chrom(r, L_b, L_s, chr)
      } else if (type == "permute") {
        unseg_permute_within_chrom(r, L_b, L_s, chr)
      }
      r_prime
    })
    y_prime <- do.call(c, r_primes)
  } else {
    # TODO: code assumes sorted 'y'
    stopifnot(all(y == GenomicRanges::sort(y)))
    # TODO: do we have to worry about missing seqlevels?
    # L_s is now a vector of the chromosome lengths
    L_s <- chrom_lens
    y_prime <- if (type == "bootstrap") {
      unseg_bootstrap_across_chrom(y, L_b, L_s)
    } else if (type == "permute") {
      unseg_permute_across_chrom(y, L_b, L_s)
    }
  }
  y_prime
}

# Block bootstrap GRanges within chromosome
#
# @param y the input GRanges
# @param L_b the length of the blocks
# @param L_s the length of the segment (chromosome)
# @param chr the name of the chromosome
unseg_bootstrap_within_chrom <- function(y, L_b, L_s=NULL, chr=NULL) {
  if (is.null(chr)) {
    chr <- as.character(seqnames(y)[1])
  }
  if (is.null(L_s)) {
    L_s <- seqlengths(y)[[chr]]
  }
  # blocks allowed to go over L_s
  n <- ceiling(L_s / L_b)
  # these blocks are the 'bait' to capture features in 'y'
  random_blocks <- IRanges::IRanges(start = round(runif(n, 1, (n - 1) * L_b + 1)),
                                    width = L_b)
  # where those blocks will move to
  rearranged_blocks <- IRanges::successiveIRanges(width(random_blocks))
  # the shift needed to move them
  block_shift <- start(rearranged_blocks) - start(random_blocks)
  # use the bait to sample features in 'y'
  fo <- IRanges::findOverlaps(random_blocks, ranges(y))
  # shift the ranges in those bait blocks
  # NOTE: we use suppress warnings, we will immediately perform a trim() operation
  suppressWarnings({
    y_prime <- shift(y[subjectHits(fo)], block_shift[queryHits(fo)])
  })
  y_prime <- trim(y_prime)
  y_prime
}

# Permute blocks of GRanges within chromosome
#
# @param y the input GRanges
# @param L_b the length of the blocks
# @param L_s the length of the segment (chromosome)
# @param chr the name of the chromosome
unseg_permute_within_chrom <- function(y, L_b, L_s=NULL, chr=NULL) {
  if (is.null(chr)) {
    chr <- as.character(seqnames(y)[1])
  }
  if (is.null(L_s)) {
    L_s <- seqlengths(y)[[chr]]
  }
  # blocks allowed to go over L_s
  n <- ceiling(L_s / L_b)
  # blocks tiling the chromosome
  blocks <- IRanges::successiveIRanges(rep(L_b, n))
  # which block do features in 'y' fall into?
  mcols(y)$block <- IRanges::findOverlaps(ranges(y), blocks, select = "first")
  # permutation order
  perm <- sample(n)
  # where those blocks will move to
  permuted_blocks <- IRanges::successiveIRanges(width(blocks))[perm]
  # the shift needed to move them
  block_shift <- start(permuted_blocks) - start(blocks)
  # shift the features in 'y'
  # NOTE: we use suppress warnings, we will immediately perform a trim() operation  
  suppressWarnings({
    y_prime <- shift(y, block_shift[mcols(y)$block])
  })
  y_prime <- trim(y_prime)
  y_prime
}

# Block bootstrap GRanges across chromosome
#
# @param y the input GRanges
# @param L_b the length of the blocks
# @param L_s the lengths of the chromosomes
unseg_bootstrap_across_chrom <- function(y, L_b, L_s) {
  # blocks allowed to go over L_s
  n_per_chrom  <- ceiling(L_s / L_b)
  # total number of tiling blocks
  n <- sum(n_per_chrom)
  # sample according to chromosome length
  random_chr <- sample(names(L_s), n, replace=TRUE, prob=L_s)
  # random starts within the tiled regions
  random_start <- round(runif(n, 1, (n_per_chrom[random_chr] - 1) * L_b + 1))
  # these blocks are the 'bait' for capturing features in 'y'
  random_blocks <- GenomicRanges::GRanges(seqnames=random_chr,
                           ranges=IRanges::IRanges(start = random_start, width = L_b))
  # where those blocks will move to
  rearranged_blocks <- GenomicRanges::tileGenome(seqlengths=L_s,
                                                 tilewidth=L_b,
                                                 cut.last.tile.in.chrom=TRUE)
  # use the bait to sample features in 'y'
  fo <- GenomicRanges::findOverlaps(random_blocks, y)
  # y has been sampled multiple times
  y_mult_hits <- y[subjectHits(fo)]
  # label which 'bait' block each feature hit in the re-sampling
  mcols(y_mult_hits)$block <- queryHits(fo)
  # shift the ranges in those bait blocks
  shift_and_swap_chrom(y_mult_hits, seqnames(rearranged_blocks),
                       start(random_blocks), start(rearranged_blocks))
}

# Permute blocks of GRanges across chomosome
#
# @param y the input GRanges
# @param L_b the length of the blocks
# @param L_s the lengths of the chromosomes
unseg_permute_across_chrom <- function(y, L_b, L_s) {
  blocks <- GenomicRanges::tileGenome(seqlengths=L_s, tilewidth=L_b,
                       cut.last.tile.in.chrom=TRUE)
    # pass along the full seqlengths of 'y'
  seqlengths(blocks) <- seqlengths(y)
  mcols(y)$block <- GenomicRanges::findOverlaps(y, blocks, select = "first")
  perm <- sample(length(blocks))
  rearranged_blocks <- blocks[perm]
  # this operation loses some ranges:
  # those that are mapped to permuted blocks that
  # are cut by `cut.last.tile.in.chrom`
  shift_and_swap_chrom(y, seqnames(rearranged_blocks),
                       start(blocks), start(rearranged_blocks))
}

# function moves featues in 'y' that fall into 'blocks'
# to the new block locations in 'rearranged_blocks', e.g.:
#
# 'y' in blocks --> 'y_prime' in rearranged_blocks
#
# and will also change the seqnames to 'chr_names'
shift_and_swap_chrom <- function(y, chr_names,
                                 random_blocks_start,
                                 rearranged_blocks_start) {
  stopifnot(length(rearranged_blocks_start) == length(random_blocks_start))
  block_shift <- rearranged_blocks_start - random_blocks_start
  idx <- mcols(y)$block
  chr_prime <- chr_names[idx]
  # this temporarily creates out-of-bound ranges
  # NOTE: we use suppress warnings, we will immediately perform a trim() operation  
  suppressWarnings({
    y_prime <- shift(y, block_shift[idx])
    seqnames(y_prime) <- chr_prime
  })
  y_prime <- trim(y_prime)
  y_prime
}
