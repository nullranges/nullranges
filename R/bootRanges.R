#' Block bootstrap genomic ranges
#'
#' Performs a block bootstrap (R times) optionally with respect
#' to a genome segmentation. Returns a bootRanges object, which
#' is essentially a GRangesList of length R.
#'
#' @param x the input GRanges
#' @param seg the segmentation GRanges, with a column ("state")
#' indicating segmentation state (optional)
#' @param blockLength the length of the blocks (for proportional blocks, this
#' is the maximal length of block)
#' @param R the number of bootstrap samples to generate
#' @param exclude the GRanges of excluded regions (optional)
#' @param excludeOption whether to "drop" or "trim" bootstrap
#' ranges that overlap a exclude region
#' @param proportionLength for segmented block bootstrap,
#' whether to use scaled block length, (scaling by the proportion
#' of the segmentation state out of the total genome length)
#' @param type the type of null generation (unsegmented only)
#' @param withinChrom whether to re-sample (bootstrap) ranges
#' across chromosomes (default) or only within chromosomes (unsegmented only)
#'
#' @return bootRanges: a GRangesList of length R with the bootstrapped ranges
#'
#' @importFrom IRanges overlapsAny
#' @importFrom stats as.formula binomial kmeans predict
#' quantile rbinom rnorm runif terms
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges successiveIRanges mid
#' @importFrom GenomicRanges tileGenome sort GRangesList
#' @importFrom GenomeInfoDb seqlengths seqlengths<- seqlevels sortSeqlevels
#'
#' @examples
#'
#' set.seed(1)
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges(0:4 * 10 + 1, width=5),
#'               seqlengths=c(chr1=50))
#' br <- bootRanges(gr, blockLength=10)
#' 
#' @export
bootRanges <- function(x, seg = NULL, blockLength, R = 1,
                       exclude = NULL,
                       excludeOption = c("drop", "trim"),
                       proportionLength = TRUE,
                       type = c("bootstrap", "permute"),
                       withinChrom = FALSE) {
  chr_lens <- seqlengths(x)
  stopifnot(all(!is.na(chr_lens)))
  excludeOption <- match.arg(excludeOption)
  if (excludeOption == "trim") {
    stopifnot(all(strand(exclude) == "*"))
  }

  br <- replicate(R, {
    if (!is.null(seg)) {
      x_prime <- seg_bootstrap(
        x = x, seg = ranges(seg),
        state = seg$state, L_b = blockLength,
        chr_names = seqnames(seg),
        chr_lens = chr_lens,
        proportion_length = proportionLength
      )
    } else {
      x_prime <- unseg_bootstrap(
        x = x, L_b = blockLength, type = type,
        within_chrom = withinChrom
      )
    }

    # deal with exclude regions:
    if (!is.null(exclude)) {
      if (excludeOption == "drop") {
        x_prime <- x_prime[!overlapsAny(x_prime, exclude, type = "any")]
      } else if (excludeOption == "trim") {
        # TODO: place outside of this function, doing gaps once?
        gap <- gaps(exclude, end = seqlengths(x))
        gap <- plyranges::filter(gap, strand == "*")
        x_prime <- plyranges::join_overlap_intersect(x_prime, gap)
      }
    }

    # sort outgoing ranges?
    GenomicRanges::sort(x_prime)
  })


  br <- GRangesList(br)
  mcols(br)$blockLength <- Rle(as.integer(blockLength))
  mcols(br)$iter <- seq_len(R)
  new("bootRanges", br)
}

# Segmentation block bootstrap sub-function
#
# @param x the input GRanges
# @param seg the segmentation ranges
# @param state the segmentation state (along seg)
# @param L_b the length of the blocks
# @param chr_names the chromosomes of the segmentation ranges
# @param chr_lens the chromosome lengths
# @param proportion_length whether to use scaled block length
#        (scaling by the proportion of the segmentation state
#         out of the total genome length)
seg_bootstrap <- function(x, seg, state, L_b,
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
    levels = seqlevels(x)
  )

  # use the bait to sample features in 'x'
  fo <- findOverlaps(random_blocks, x)
  # x has been sampled multiple times
  x_mult_hits <- x[subjectHits(fo)]
  # label which 'bait' block each feature hit in the re-sampling
  mcols(x_mult_hits)$block <- queryHits(fo)
  # shift the ranges in those bait blocks
  shift_and_swap_chrom(
    x_mult_hits, rearr_blocks_chr_names,
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
