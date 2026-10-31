# bootRanges() as nullranges shipped it before its bulk and compiled paths
# (R/bootRanges.R and R/bootUnsegmented.R of 1.19.2), for the inputs where
# those paths could differ from it: bootRanges() runs it for a 'y' that
# fails GRanges' validity, an RNG whose state .Random.seed does not hold, and
# any call on which its own code signals a warning or an error, after putting
# back .Random.seed as it was before the first draw. The text below is that
# code unchanged, with the roxygen block removed; it is defined inside local()
# so that its functions keep their names, and the calls their errors and
# warnings carry, without masking the package's own.
bootRanges_base <- local({

bootRanges <- function(y, blockLength, R = 1,
                       seg = NULL,
                       exclude = NULL,
                       excludeOption = c("drop", "trim"),
                       proportionLength = TRUE,
                       type = c("bootstrap", "permute"),
                       withinChrom = FALSE,
                       storeBlockLength = FALSE) {
  chr_lens <- seqlengths(y)
  stopifnot(all(!is.na(chr_lens)))
  excludeOption <- match.arg(excludeOption)
  if (excludeOption == "trim") {
    stopifnot(all(strand(exclude) == "*"))
  }
  if (!is.null(seg)) {
    stopifnot("state" %in% names(mcols(seg)))
  }

  br <- replicate(R, {
    if (!is.null(seg)) {
      y_prime <- seg_bootstrap(
        y = y, seg = ranges(seg),
        state = seg$state, L_b = blockLength,
        chr_names = seqnames(seg),
        chr_lens = chr_lens,
        proportion_length = proportionLength
      )
    } else {
      y_prime <- unseg_bootstrap(
        y = y, L_b = blockLength, type = type,
        within_chrom = withinChrom
      )
    }

    # deal with exclude regions:
    if (!is.null(exclude)) {
      if (excludeOption == "drop") {
        y_prime <- y_prime[!overlapsAny(y_prime, exclude, type = "any")]
      } else if (excludeOption == "trim") {
        # TODO: place outside of this function, doing gaps once?
        gap <- gaps(exclude, end = seqlengths(y))
        gap <- gap[strand(gap) == "*"]
        y_prime <- plyranges::join_overlap_intersect(y_prime, gap)
      }
    }
    y_prime
  })

  lens <- lengths(br)
  br <- do.call(c, br)
  mcols(br)$iter <- Rle(factor(rep(seq_len(R), lens), levels=seq_len(R)))
  if (storeBlockLength) {
    mcols(br)$blockLength <- Rle(as.integer(blockLength))
  }
  mcols(br)$block <- NULL # this is not needed except for debugging perhaps
  new("BootRanges", br)
}

# Segmentation block bootstrap sub-function
#
# @param y the input GRanges
# @param seg the segmentation ranges
# @param state the segmentation state (along seg)
# @param L_b the length of the blocks
# @param chr_names the chromosomes of the segmentation ranges
# @param chr_lens the chromosome lengths
# @param proportion_length whether to use scaled block length
#        (scaling by the proportion of the segmentation state
#         out of the total genome length)
seg_bootstrap <- function(y, seg, state, L_b,
                          chr_names, chr_lens,
                          proportion_length = TRUE) {

  # sequence along the states: 1,2,..,S
  esses <- seq_len(max(state))
  chr_names <- as.character(chr_names)

  # break apart this into two sub-functions, defined below
  if (proportion_length) {
    out <- seg_bootstrap_prop(seg, state, L_b, chr_names, chr_lens, esses)
  } else {
    out <- seg_bootstrap_noprop(seg, state, L_b, chr_names, chr_lens, esses)
  }

  random_blocks <- out$random_blocks
  rearr_blocks_chr_names <- out$rearr_blocks_chr_names
  random_blocks_start <- out$random_blocks_start
  rearranged_blocks_start <- out$rearranged_blocks_start

  rearr_blocks_chr_names <- factor(rearr_blocks_chr_names,
    levels = seqlevels(y)
  )

  # use the bait to sample features in 'y'
  fo <- findOverlaps(random_blocks, y)
  # y has been sampled multiple times
  y_mult_hits <- y[subjectHits(fo)]
  # label which 'bait' block each feature hit in the re-sampling
  mcols(y_mult_hits)$block <- queryHits(fo)
  # shift the ranges in those bait blocks
  shift_and_swap_chrom(
    y_mult_hits, rearr_blocks_chr_names,
    random_blocks_start, rearranged_blocks_start
  )
}

# breaking up the seg_bootstrap() code... proportional length case
seg_bootstrap_prop <- function(seg, state, L_b, chr_names, chr_lens, esses) {
  # length of genome
  L_g <- sum(chr_lens)
  # segmentation state total length
  L_s <- unlist(lapply(esses, function(s) {
    sum(width(seg)[state == s])
  }), use.names = FALSE)
  # the block width for each segmentation state is scaled
  # down to the segmentation state size, e.g. if segmentation state
  # is half of the genome, then the block width is half of L_b
  L_b_per_state <- round(L_b * L_s / L_g)

  # loops over states 's'
  obj <- lapply(esses, function(s) {
    # the length of the blocks for this state
    L_bs <- L_b_per_state[s]
    # the segmentation ranges for this state
    seg_s <- seg[state == s]
    # number of blocks within each range of the segmentation
    nblocks_per_rng <- ceiling(width(seg_s) / L_bs)
    # total number of tiling blocks
    n <- sum(nblocks_per_rng)
    # sample ranges according to their length
    random_rng <- sample(length(seg_s), n,
      replace = TRUE,
      prob = width(seg_s)
    )
    # random starts within the ranges, as spec by 'nblocks per rng'
    rmin <- start(seg_s)[random_rng]
    rmax <- (nblocks_per_rng[random_rng] - 1) * L_bs + rmin
    random_start_s <- round(runif(n, rmin, rmax))
    # names of chromosomes for random ranges of state 's'
    random_chr_names_s <- chr_names[state == s][random_rng]
    # start positions of blocks tiling ranges of state 's'
    rearr_blocks_start_s <- unlist(lapply(seq_along(seg_s), function(j) {
      seq(from = start(seg_s)[j], to = end(seg_s)[j], by = L_bs)
    }), use.names = FALSE)
    # names of chromosomes for blocks tiling ranges of state 's'
    rearr_blocks_chr_names_s <- rep(
      chr_names[state == s], nblocks_per_rng
    )
    return(list(
      random_start_s, random_chr_names_s,
      rearr_blocks_start_s, rearr_blocks_chr_names_s
    ))
  })

  # number of blocks for each segmentation state
  nblocks_per_state <- lengths(lapply(obj, `[[`, 1))

  # unlist the results in 'obj'
  random_blocks_start <- unlist(lapply(obj, `[[`, 1))
  random_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 2))
  rearranged_blocks_start <- unlist(lapply(obj, `[[`, 3))
  rearr_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 4))

  # construct the random blocks
  random_blocks <- GRanges(
    seqnames = random_blocks_chr_names,
    ranges = IRanges(
      start = random_blocks_start,
      width = rep(
        L_b_per_state,
        nblocks_per_state
      )
    )
  )

  list(
    random_blocks = random_blocks,
    rearr_blocks_chr_names = rearr_blocks_chr_names,
    random_blocks_start = random_blocks_start,
    rearranged_blocks_start = rearranged_blocks_start
  )
}

# breaking up the seg_bootstrap() code... non-proportional length case
seg_bootstrap_noprop <- function(seg, state, L_b, chr_names, chr_lens, esses) {

  # number of rearranged blocks within each range of the segmentation
  nblocks_per_seg <- ceiling(width(seg) / L_b)
  # total number of tiling blocks
  n <- sum(nblocks_per_seg)
  # sample segments according to segmentation length
  random_seg <- sample(length(seg), n, replace = TRUE, prob = width(seg))
  # random starts within segments, as specified by 'nblocks_per_seg'
  rmin <- start(seg)[random_seg]
  rmax <- (nblocks_per_seg[random_seg] - 1) * L_b + rmin
  random_start <- round(runif(n, rmin, rmax))
  # start positions of rearranged blocks tiling the segments
  rearr_blocks_start_list <- lapply(seq_along(seg), function(j) {
    seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
  })

  # loop over states 's', splitting out the random segments for 's'
  obj <- lapply(esses, function(s) {
    # index of those random segments chosen that are state 's'
    idx <- which(state[random_seg] == s)
    # random starts for state s
    random_start_s <- random_start[idx]
    # names of chromosomes for random segments that are state 's'
    random_chr_names_s <- chr_names[random_seg[idx]]
    # start positions of blocks tiling segment for state 's'
    rearr_blocks_start_s <- unlist(
      rearr_blocks_start_list[state == s],
      use.names = FALSE
    )
    # names of chromosomes for blocks tilings segment for state 's'
    rearr_blocks_chr_names_s <- rep(
      chr_names[state == s],
      nblocks_per_seg[state == s]
    )
    # see if we have too few or too many blocks for this state
    rand_n <- length(random_start_s)
    rearr_n <- length(rearr_blocks_start_s)
    # if we have too few random blocks, up-sample them
    if (rand_n < rearr_n) {
      rand_idx <- sample(rand_n, rearr_n, replace = TRUE)
    } else {
      # otherwise down-sample to the number of re-arranged blocks
      # (or stay if the number happened to be equal)
      rand_idx <- seq_len(rearr_n)
    }
    return(list(
      random_start_s[rand_idx],
      random_chr_names_s[rand_idx],
      rearr_blocks_start_s,
      rearr_blocks_chr_names_s
    ))
  })

  # unlist the results in 'obj'
  random_blocks_start <- unlist(lapply(obj, `[[`, 1))
  random_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 2))
  rearranged_blocks_start <- unlist(lapply(obj, `[[`, 3))
  rearr_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 4))

  # construct the random blocks
  random_blocks <- GRanges(
    seqnames = random_blocks_chr_names,
    ranges = IRanges(
      start = random_blocks_start,
      width = L_b
    )
  )

  list(
    random_blocks = random_blocks,
    rearr_blocks_chr_names = rearr_blocks_chr_names,
    random_blocks_start = random_blocks_start,
    rearranged_blocks_start = rearranged_blocks_start
  )
}

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
    # TODO: what happens when no features occur on a chromosome?
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
  shift_and_swap_chrom(
    y_mult_hits, seqnames(rearranged_blocks),
    start(random_blocks), start(rearranged_blocks)
  )
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
  # this temporarily creates out-of-bound ranges
  # NOTE: we use suppress warnings, we will immediately perform a trim() operation  
  suppressWarnings({
    y_prime <- shift(y, block_shift[idx])
    seqnames(y_prime) <- chr_prime
  })
  y_prime <- trim(y_prime)
  y_prime <- y_prime[width(y_prime) > 0]
  y_prime
}

bootRanges
})
