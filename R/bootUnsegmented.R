# The part of the un-segmented bootstrap that does not depend on the
# random draws: checks, chromosome table and the tiling blocks.
# bootRanges() calls it once and passes the result to every replicate.
unseg_bootstrap_setup <- function(y, L_b,
                                  type = c("bootstrap", "permute"),
                                  within_chrom = FALSE) {
  type <- match.arg(type)
  chrom_lens <- seqlengths(y)
  # TODO: what kind of solution do we want for this?
  # it's probably not a good idea to have L_b larger than one of the sequences
  stopifnot(all(chrom_lens >= L_b))
  tab <- table(seqnames(y))
  chroms <- names(tab)
  blocks <- NULL
  if (!within_chrom) {
    # TODO: code assumes sorted 'y'
    # is.unsorted() tests the keys sort() orders a GRanges by, so when it
    # is FALSE sort(y) equals y and the check passes
    if (!identical(as.vector(class(y)), "GRanges") || is.unsorted(y)) {
      stopifnot(all(y == GenomicRanges::sort(y)))
    }
    # TODO: what happens when no features occur on a chromosome?
    blocks <- if (type == "bootstrap") tile_starts_genome(chrom_lens, L_b)
    if (is.null(blocks)) {
      blocks <- GenomicRanges::tileGenome(seqlengths=chrom_lens,
                                          tilewidth=L_b,
                                          cut.last.tile.in.chrom=TRUE)
    }
    if (type == "permute") {
      # pass along the full seqlengths of 'y'
      seqlengths(blocks) <- seqlengths(y)
      mcols(y)$block <- GenomicRanges::findOverlaps(y, blocks, select = "first")
    }
  }
  list(y = y, L_b = L_b, type = type, within_chrom = within_chrom,
       chrom_lens = chrom_lens, chroms = chroms, blocks = blocks)
}

# start() and seqnames() of tileGenome(seqlengths = chrom_lens, tilewidth =
# L_b, cut.last.tile.in.chrom = TRUE), all the bootstrap reads of it, as a
# list: with a whole L_b a sequence of length len is tiled from 1 + L_b * k,
# k = 0 .. ceiling(len / L_b) - 1. NULL wherever tileGenome() could fail or
# cut otherwise: L_b not a single whole number >= 1, no sequence, names
# missing, empty or duplicated, a length missing or below L_b.
#
# @param chrom_lens the chromosome lengths, named
# @param L_b the length of the blocks
tile_starts_genome <- function(chrom_lens, L_b) {
  nm <- names(chrom_lens)
  if (!is.numeric(L_b) || length(L_b) != 1L || is.object(L_b) ||
      !is.finite(L_b) || L_b < 1 || L_b != floor(L_b) ||
      !is.integer(chrom_lens) || length(chrom_lens) == 0L || is.null(nm) ||
      anyNA(nm) || any(nm == "") || anyDuplicated(nm) ||
      anyNA(chrom_lens) || any(chrom_lens < L_b)) {
    return(NULL)
  }
  n_tile <- as.integer(ceiling(unname(chrom_lens) / L_b))
  list(
    start = as.integer(1 + L_b * (sequence(n_tile) - 1L)),
    seqnames = Rle(factor(nm, levels = nm), n_tile)
  )
}

# start() and seqnames() of the tiling blocks: tileGenome()'s GRanges, or
# the list tile_starts_genome() gives in its place
#
# @param blocks the blocks tiling the genome
tile_parts <- function(blocks) {
  if (is.list(blocks)) {
    return(blocks)
  }
  list(start = start(blocks), seqnames = seqnames(blocks))
}

# One un-segmented bootstrap (or permutation) replicate
#
# @param setup the output of unseg_bootstrap_setup()
unseg_bootstrap <- function(setup) {
  y <- setup$y
  L_b <- setup$L_b
  type <- setup$type
  chrom_lens <- setup$chrom_lens
  # first, simple case: sampling or permuting within chrom

  # TODO this loses the metadata columns right now...

  if (setup$within_chrom) {
    r_primes <- lapply(setup$chroms, function(chr) {
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
    # L_s is now a vector of the chromosome lengths
    L_s <- chrom_lens
    y_prime <- if (type == "bootstrap") {
      unseg_bootstrap_across_chrom(y, L_b, L_s, setup$blocks)
    } else if (type == "permute") {
      unseg_permute_across_chrom(y, setup$blocks)
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
# @param rearranged_blocks the blocks tiling the genome (tileGenome),
# where the bait blocks will move to
unseg_bootstrap_across_chrom <- function(y, L_b, L_s, rearranged_blocks) {
  blocks <- unseg_across_draws(1L, L_b, L_s, rearranged_blocks)
  move_blocks(y, blocks, blocks$draws[[1L]])
}

# The random bait blocks of 'n_rep' replicates of the bootstrap across
# chromosomes, in replicate order, and the blocks they move to; the
# value is laid out as seg_bootstrap_draws()'s
#
# @param n_rep the number of replicates
# @param L_b the length of the blocks
# @param L_s the lengths of the chromosomes
# @param rearranged_blocks the blocks tiling the genome (tileGenome),
# where the bait blocks will move to
# @param compiled whether seg_draws_compiled() may make the draws, as one
# table over the chromosomes; then 'draws' is NULL and 'bait_start',
# 'bait_seg' (the chromosome's index in L_s), 'bait_len' and 'chr_names'
# hold them, as seg_prop_compiled() lays them out
unseg_across_draws <- function(n_rep, L_b, L_s, rearranged_blocks,
                               compiled = FALSE) {
  # blocks allowed to go over L_s
  n_per_chrom  <- ceiling(L_s / L_b)
  # total number of tiling blocks
  n <- sum(n_per_chrom)
  # the loop's draws per replicate: sample() over the chromosomes, then
  # runif() up to (n_per_chrom - 1) * L_b + 1
  d <- if (compiled) {
    seg_draws_compiled(list(list(
      index = seq_along(L_s), width = L_s, start = rep(1L, length(L_s)),
      nblocks = n_per_chrom, L = L_b, n = n
    )), n_rep)
  }
  tiles <- tile_parts(rearranged_blocks)
  if (!is.null(d)) {
    return(list(
      draws = NULL, width = L_b,
      target_start = tiles$start,
      target_chr = tiles$seqnames,
      bait_start = d[[1L]], bait_seg = d[[2L]], bait_len = d[[3L]],
      chr_names = names(L_s)
    ))
  }
  draws <- lapply(seq_len(n_rep), function(r) {
    # sample according to chromosome length
    random_chr <- sample(names(L_s), n, replace=TRUE, prob=L_s)
    # random starts within the tiled regions
    random_start <- round(runif(n, 1, (n_per_chrom[random_chr] - 1) * L_b + 1))
    list(start = random_start, chr = random_chr)
  })
  list(draws = draws, width = L_b,
       target_start = tiles$start,
       target_chr = tiles$seqnames)
}

# Permute blocks of GRanges across chomosome
#
# @param y the input GRanges, with the block each range falls into
# in mcols(y)$block
# @param blocks the blocks tiling the genome (tileGenome)
unseg_permute_across_chrom <- function(y, blocks) {
  perm <- sample(length(blocks))
  rearranged_blocks <- blocks[perm]
  # this operation loses some ranges:
  # those that are mapped to permuted blocks that
  # are cut by `cut.last.tile.in.chrom`
  shift_and_swap_chrom(
    y, seqnames(rearranged_blocks),
    start(blocks), start(rearranged_blocks)
  )
}

# function moves featues in 'y' that fall into 'blocks'
# to the new block locations in 'rearranged_blocks', e.g.:
#
# 'y' in blocks --> 'y_prime' in rearranged_blocks
#
# and will also change the seqnames to 'chr_names'
#
# some features after trimming will have width 0 and are removed
shift_and_swap_chrom <- function(y, chr_names,
                                 random_blocks_start,
                                 rearranged_blocks_start) {
  stopifnot(length(rearranged_blocks_start) == length(random_blocks_start))
  block_shift <- rearranged_blocks_start - random_blocks_start
  idx <- mcols(y)$block
  chr_prime <- chr_names[idx]
  y_prime <- shift_swap_trim(y, block_shift[idx], chr_prime)
  if (is.null(y_prime)) {
    # this temporarily creates out-of-bound ranges
    # NOTE: we use suppress warnings, we will immediately perform a trim() operation
    suppressWarnings({
      y_prime <- shift(y, block_shift[idx])
      seqnames(y_prime) <- chr_prime
    })
    y_prime <- trim(y_prime)
  }
  y_prime <- y_prime[width(y_prime) > 0]
  y_prime
}

# shift(y, shift), then seqnames(y) <- chr_prime, then trim(), computed on
# the start and width vectors and written in one update without validity
# passes. Returns NULL for any input where one of those three calls could
# warn, fail or treat the object differently (not a plain GRanges, NA or
# overflowing positions, seqnames not a factor on seqlevels(y)); the caller
# then takes the original path.
shift_swap_trim <- function(y, shift, chr_prime) {
  n <- length(y)
  if (n == 0L || !identical(as.vector(class(y)), "GRanges") ||
      !is.numeric(shift) || length(shift) != n) {
    return(NULL)
  }
  shift <- suppressWarnings(as.integer(shift))
  if (anyNA(shift)) {
    return(NULL)
  }
  new_seqnames <- if (is(chr_prime, "Rle")) chr_prime else Rle(chr_prime)
  if (length(new_seqnames) != n || !is.factor(runValue(new_seqnames)) ||
      anyNA(runValue(new_seqnames)) ||
      !identical(levels(new_seqnames), seqlevels(y))) {
    return(NULL)
  }
  rng <- ranges(y)
  new_width <- width(rng)
  new_start <- suppressWarnings(start(rng) + shift)
  new_end <- suppressWarnings(new_start + (new_width - 1L))
  if (anyNA(new_start) || anyNA(new_end)) {
    return(NULL)
  }
  # trim(): clip to [1, seqlength] unless the seqlevel is circular or has
  # no length, as get_out_of_bound_index() and restrict() do
  seqlen <- unname(seqlengths(y))
  bounded <- !(unname(Seqinfo::isCircular(y)) %in% TRUE | is.na(seqlen))
  code <- as.integer(new_seqnames)
  chr_len <- seqlen[code]
  oob <- which(bounded[code] & (new_start < 1L | new_end > chr_len))
  if (length(oob) != 0L) {
    oob_start <- pmax(new_start[oob], 1L)
    oob_end <- pmin(new_end[oob], chr_len[oob])
    new_start[oob] <- oob_start
    new_width[oob] <- pmax(oob_end - oob_start + 1L, 0L)
  }
  new_ranges <- update_ranges(rng, start = new_start, width = new_width)
  update(y, ranges = new_ranges, seqnames = new_seqnames, check = FALSE)
}
