#' Block bootstrap for genomic ranges
#'
#' Performs a block bootstrap, optionally with respect
#' to a genome segmentation. Returns a \code{bootRanges} object,
#' which is a GRanges object with all the ranges concatenated,
#' and iteration and block length indicated by metadata columns
#'
#' @param y the GRanges to bootstrap sample
#' @param blockLength the length of the blocks
#' (for proportional blocks, this is the maximal length of a block)
#' @param R the number of bootstrap samples to generate
#' @param seg the segmentation GRanges, with a column ("state")
#' indicating segmentation state (optional)
#' @param exclude the GRanges of excluded regions (optional)
#' @param excludeOption whether to \code{"drop"} or \code{"trim"}
#' bootstrap ranges that overlap a excluded region
#' @param proportionLength for the segmented block bootstrap,
#' whether to use scaled block lengths, (scaling by the proportion
#' of the segmentation state out of the total genome length)
#' @param type the type of null generation (un-segmented bootstrap only)
#' @param withinChrom whether to re-sample (bootstrap) ranges
#' across chromosomes (default) or only within chromosomes
#' (un-segmented bootstrap only)
#'
#' @return a BootRanges (GRanges object) with the bootstrapped ranges,
#' where iteration and block length are recorded as metadata columns
#'
#' @importFrom IRanges overlapsAny
#' @importFrom stats as.formula binomial kmeans predict
#' quantile rbinom rnorm runif terms
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges successiveIRanges mid
#' @importFrom GenomicRanges tileGenome sort GRangesList
#' @importFrom GenomeInfoDb seqlengths seqlengths<- seqlevels sortSeqlevels
#'
#' @references
#'
#' bootRanges manuscript:
#'
#' Wancen Mu, Eric S. Davis, Stuart Lee, Mikhail G. Dozmorov,
#' Douglas H. Phanstiel, Michael I. Love.
#' 2022. "bootRanges: Flexible generation of null sets
#' of genomic ranges for hypothesis testing."
#' bioRxiv. doi: 10.1101/2022.09.02.506382
#' 
#' Original method describing the segmented block bootstrap for genomic features:
#' 
#' Bickel, Peter J., Nathan Boley, James B. Brown,
#' Haiyan Huang, and Nancy R. Zhang.
#' 2010. "Subsampling Methods for Genomic Inference."
#' The Annals of Applied Statistics 4 (4): 1660–97.
#' doi: 10.1214/10-AOAS363
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
bootRanges <- function(y, blockLength, R = 1,
                       seg = NULL,
                       exclude = NULL,
                       excludeOption = c("drop", "trim"),
                       proportionLength = TRUE,
                       type = c("bootstrap", "permute"),
                       withinChrom = FALSE) {
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
        gap <- plyranges::filter(gap, strand == "*")
        y_prime <- plyranges::join_overlap_intersect(y_prime, gap)
      }
    }
    y_prime
  })

  lens <- lengths(br)
  br <- do.call(c, br)
  mcols(br)$iter <- Rle(factor(rep(seq_len(R), lens), levels=seq_len(R)))
  mcols(br)$blockLength <- Rle(as.integer(blockLength))
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
