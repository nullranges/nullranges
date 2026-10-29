#' Block bootstrap for genomic ranges
#'
#' Performs a block bootstrap, optionally with respect
#' to a genome segmentation. Returns a `bootRanges` object,
#' which is a GRanges object with all the ranges from the
#' bootstrap iterations concatenated.
#' 
#' Note that this function requires input ranges have associated
#' `seqlengths`, and that these must not be shorter than `blockLength`.
#' See `Seqinfo`, `seqlevels`, and `keepStandardChromosomes` functions
#' and their use in the Quick Start section of the vignette.
#'
#' @param y the GRanges to bootstrap sample
#' @param blockLength the length of the blocks
#' (for proportional blocks, this is the maximal length of a block)
#' @param R the number of bootstrap samples to generate
#' @param seg the segmentation GRanges, with a column ("state")
#' indicating segmentation state (optional)
#' @param exclude the GRanges of excluded regions (optional)
#' @param excludeOption whether to `"drop"` or `"trim"`
#' bootstrap ranges that overlap a excluded region
#' @param proportionLength for the segmented block bootstrap,
#' whether to use scaled block lengths (scaling by the proportion
#' of the segmentation state out of the total genome length).
#' That is, the resulting blocks will be of size less than
#' `blockLength`
#' @param type the type of null generation (un-segmented bootstrap only)
#' @param withinChrom whether to re-sample (bootstrap) ranges
#' across chromosomes (default) or only within chromosomes
#' (un-segmented bootstrap only)
#' @param storeBlockLength whether to save blockLength as a
#' metadata column
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
#' @importFrom Seqinfo seqlengths seqlengths<- seqlevels sortSeqlevels genome
#'
#' @references
#'
#' bootRanges manuscript:
#'
#' Wancen Mu, Eric S. Davis, Stuart Lee, Mikhail G. Dozmorov,
#' Douglas H. Phanstiel, Michael I. Love.
#' 2023. "bootRanges: Flexible generation of null sets
#' of genomic ranges for hypothesis testing."
#' Bioinformatics. doi: 10.1093/bioinformatics/btad190
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

  # done once for all replicates; skipped, with its checks, when R < 1
  if (is.null(seg) && R >= 1) {
    unseg_setup <- unseg_bootstrap_setup(
      y = y, L_b = blockLength, type = type,
      within_chrom = withinChrom
    )
  }

  # as many replicates as replicate(R, ...) makes
  n_rep <- length(integer(R))
  bait <- n_rep >= 1L && (!is.null(seg) ||
    (unseg_setup$type == "bootstrap" && isFALSE(unseg_setup$within_chrom)))

  if (bait) {
    # the random bait blocks of every replicate first, in replicate order,
    # then the features they capture are moved for all replicates at once
    blocks <- if (!is.null(seg)) {
      seg_bootstrap_draws(
        n_rep, y = y, seg = ranges(seg),
        state = seg$state, L_b = blockLength,
        chr_names = seqnames(seg),
        chr_lens = chr_lens,
        proportion_length = proportionLength
      )
    } else {
      unseg_across_draws(
        n_rep, unseg_setup$L_b, unseg_setup$chrom_lens, unseg_setup$blocks
      )
    }
    # NULL for any input where the replicates, run one at a time, could
    # warn, fail or differ: then they are run one at a time
    bulk <- tryCatch(
      move_blocks_bulk(y, blocks, exclude, excludeOption),
      warning = function(w) NULL, error = function(e) NULL
    )
    if (!is.null(bulk)) {
      lens <- bulk$lens
      br <- bulk$ranges
    } else {
      br <- lapply(blocks$draws, function(d) {
        exclude_ranges(move_blocks(y, blocks, d), y, exclude, excludeOption)
      })
      lens <- lengths(br)
      br <- do.call(c, br)
    }
  } else {
    br <- replicate(R, {
      y_prime <- unseg_bootstrap(unseg_setup)
      exclude_ranges(y_prime, y, exclude, excludeOption)
    })
    lens <- lengths(br)
    br <- do.call(c, br)
  }

  mcols(br)$iter <- Rle(factor(rep(seq_len(R), lens), levels=seq_len(R)))
  if (storeBlockLength) {
    mcols(br)$blockLength <- Rle(as.integer(blockLength))
  }
  mcols(br)$block <- NULL # this is not needed except for debugging perhaps
  new("BootRanges", br)
}

# deal with exclude regions, for one replicate
exclude_ranges <- function(y_prime, y, exclude, excludeOption) {
  if (!is.null(exclude)) {
    if (excludeOption == "drop") {
      y_prime <- y_prime[!overlapsAny(y_prime, exclude, type = "any")]
    } else if (excludeOption == "trim") {
      gap <- gaps(exclude, end = seqlengths(y))
      gap <- gap[strand(gap) == "*"]
      y_prime <- plyranges::join_overlap_intersect(y_prime, gap)
    }
  }
  y_prime
}

# The features of 'y' under one replicate's bait blocks 'd', moved to the
# blocks' targets
#
# @param y the input GRanges
# @param blocks the output of seg_bootstrap_draws() or unseg_across_draws()
# @param d one replicate's bait blocks, an element of blocks$draws
move_blocks <- function(y, blocks, d) {
  random_blocks <- GRanges(
    seqnames = d$chr,
    ranges = IRanges(start = d$start, width = blocks$width)
  )
  # use the bait to sample features in 'y'
  fo <- findOverlaps(random_blocks, y)
  # y has been sampled multiple times
  y_mult_hits <- y[subjectHits(fo)]
  # label which 'bait' block each feature hit in the re-sampling
  mcols(y_mult_hits)$block <- queryHits(fo)
  # shift the ranges in those bait blocks
  shift_and_swap_chrom(
    y_mult_hits, blocks$target_chr,
    start(random_blocks), blocks$target_start
  )
}

# move_blocks() and exclude_ranges() for every replicate in one pass: the
# bait blocks of all replicates in one GRanges, the moved features in one
# GRanges, whose 'block' column holds the replicate. Returns the ranges
# and the number per replicate, or NULL where one replicate at a time
# could differ (a replicate with fewer bait blocks than targets, a bait
# block on a sequence 'y' lacks or past the end of a circular one, an
# input shift_swap_trim() declines). The moved positions are computed on
# vectors and 'y' is subset once, for the ranges kept.
move_blocks_bulk <- function(y, blocks, exclude, excludeOption) {
  draws <- blocks$draws
  n_rep <- length(draws)
  n_blocks <- length(blocks$target_start)
  bait_start <- lapply(draws, `[[`, "start")
  if (any(lengths(bait_start) != n_blocks)) {
    return(NULL)
  }
  bait_start <- unlist(bait_start, use.names = FALSE)
  bait_chr <- unlist(lapply(draws, `[[`, "chr"), use.names = FALSE)
  bait_code <- match(bait_chr, seqlevels(y))
  circular <- unname(Seqinfo::isCircular(y)) %in% TRUE
  if (anyNA(bait_code) || any(circular[bait_code] & (bait_start < 1 |
      bait_start > unname(seqlengths(y))[bait_code]))) {
    return(NULL)
  }
  random_blocks <- GRanges(
    seqnames = bait_chr,
    ranges = IRanges(
      start = bait_start,
      width = rep_len(blocks$width, length(bait_start))
    )
  )
  bait_rep <- rep(seq_len(n_rep), each = n_blocks)
  fo <- bait_overlaps(random_blocks, bait_code, bait_rep, n_rep, y)
  hit_block <- fo$query
  block_shift <- rep(blocks$target_start, n_rep) - start(random_blocks)
  target_chr <- blocks$target_chr
  if (is(target_chr, "Rle")) {
    target_chr <- decode(target_chr)
  }
  moved <- shift_swap_trim_positions(
    y, fo$subject, block_shift[hit_block], target_chr, n_rep, hit_block
  )
  if (is.null(moved)) {
    return(NULL)
  }
  block <- bait_rep[hit_block]
  keep <- moved$width > 0L
  if (!is.null(exclude) && excludeOption == "drop") {
    keep <- keep & !overlapsAny(
      moved_ranges(y, moved$s, moved$code, moved$start, moved$width,
                   moved$chr_attr),
      exclude, type = "any"
    )
  }
  if (!is.null(exclude) && excludeOption == "trim") {
    gap <- gaps(exclude, end = seqlengths(y))
    gap <- gap[strand(gap) == "*"]
    i <- which(keep)
    start_i <- moved$start[i]
    width_i <- moved$width[i]
    clip <- intersect_gaps(
      y, gap, moved$code[i], start_i, start_i + (width_i - 1L)
    )
    if (is.null(clip)) {
      y_prime <- moved_ranges(
        y, moved$s[i], moved$code[i], start_i, width_i,
        moved$chr_attr, block[i]
      )
      y_prime <- plyranges::join_overlap_intersect(y_prime, gap)
    } else {
      i <- i[clip$q]
      y_prime <- moved_ranges(
        y, moved$s[i], moved$code[i], clip$start, clip$width,
        moved$chr_attr, block[i], clip$seqinfo
      )
    }
  } else {
    i <- which(keep)
    y_prime <- moved_ranges(
      y, moved$s[i], moved$code[i], moved$start[i], moved$width[i],
      moved$chr_attr, block[i]
    )
  }
  list(ranges = y_prime, lens = tabulate(mcols(y_prime)$block, n_rep))
}

# What shift_swap_trim(y[s], shift, Rle(rep(target_chr, n_rep)[hit_block]))
# computes, as vectors: the subject index 's', the seqlevel code, start
# and width of each moved range, and the attributes of the seqnames
# factor. NULL wherever shift_swap_trim() declines, and for a target_chr
# carrying attributes besides levels and class.
#
# @param y the input GRanges
# @param s the index in 'y' of each moved range
# @param shift the shift of each moved range
# @param target_chr the sequence each block moves to, a factor
# @param n_rep the number of replicates
# @param hit_block the block of each moved range, among all replicates'
shift_swap_trim_positions <- function(y, s, shift, target_chr, n_rep,
                                      hit_block) {
  n <- length(s)
  if (n == 0L || !identical(as.vector(class(y)), "GRanges") ||
      !is.numeric(shift) || length(shift) != n) {
    return(NULL)
  }
  shift <- suppressWarnings(as.integer(shift))
  if (anyNA(shift)) {
    return(NULL)
  }
  chr_attr <- attributes(target_chr)
  if (!is.factor(target_chr) ||
      !setequal(names(chr_attr), c("levels", "class")) ||
      !identical(levels(target_chr), seqlevels(y))) {
    return(NULL)
  }
  code <- rep(as.integer(target_chr), n_rep)[hit_block]
  if (anyNA(code)) {
    return(NULL)
  }
  rng <- ranges(y)
  new_width <- width(rng)[s]
  new_start <- suppressWarnings(start(rng)[s] + shift)
  new_end <- suppressWarnings(new_start + (new_width - 1L))
  if (anyNA(new_start) || anyNA(new_end)) {
    return(NULL)
  }
  # trim(): clip to [1, seqlength] unless the seqlevel is circular or has
  # no length, as get_out_of_bound_index() and restrict() do
  seqlen <- unname(seqlengths(y))
  bounded <- !(unname(Seqinfo::isCircular(y)) %in% TRUE | is.na(seqlen))
  chr_len <- seqlen[code]
  oob <- new_start < 1L | new_end > chr_len
  if (!all(bounded)) {
    oob <- bounded[code] & oob
  }
  oob <- which(oob)
  if (length(oob) != 0L) {
    oob_start <- pmax(new_start[oob], 1L)
    oob_end <- pmin(new_end[oob], chr_len[oob])
    new_start[oob] <- oob_start
    new_width[oob] <- pmax(oob_end - oob_start + 1L, 0L)
  }
  list(s = s, code = code, start = new_start, width = new_width,
       chr_attr = chr_attr)
}

# Moved ranges as one subset of 'y': y's ranges, names, strand and
# metadata columns at 's', with the replicate added as column 'block' (no
# metadata columns and no names when 'block' is NULL), the given starts,
# widths and seqlevel codes, and 'seqinfo' (y's when NULL). A strand of
# one run is not subset.
#
# @param y the input GRanges
# @param s the index in 'y' of each range
# @param code the index of each range's sequence in seqlevels(y)
# @param start,width the start and width of each range
# @param chr_attr the attributes of the seqnames factor
# @param block the replicate of each range, or NULL
# @param seqinfo the seqinfo of the result, or NULL
moved_ranges <- function(y, s, code, start, width, chr_attr, block = NULL,
                         seqinfo = NULL) {
  attributes(code) <- chr_attr
  n <- length(s)
  ans_strand <- strand(y)
  one_strand <- nrun(ans_strand) == 1L
  if (!one_strand || !is.null(block)) {
    # checked once for every slot subset, as y[s] does
    s <- normalizeSingleBracketSubscript(s, y, as.NSBS = TRUE)
  }
  ans_strand <- if (one_strand) {
    Rle(runValue(ans_strand), n)
  } else {
    extractROWS(ans_strand, s)
  }
  if (is.null(block)) {
    ans_ranges <- new2("IRanges", start = start, width = width, check = FALSE)
    ans_mcols <- new2("DFrame", nrows = n, check = FALSE)
  } else {
    ans_ranges <- update_ranges(
      extractROWS(ranges(y), s), start = start, width = width
    )
    ans_mcols <- extractROWS(mcols(y, use.names = FALSE), s)
    ans_mcols$block <- block
  }
  ans <- update(y, seqnames = Rle(code), ranges = ans_ranges,
                strand = ans_strand, elementMetadata = ans_mcols,
                check = FALSE)
  if (!is.null(seqinfo)) {
    ans <- update(ans, seqinfo = seqinfo, check = FALSE)
  }
  ans
}

# plyranges::join_overlap_intersect(x, gap) for the moved ranges whose
# seqlevel codes, starts and ends are given, where 'gap' is disjoint and
# sorted along each sequence, as gaps() returns it: each range clipped to
# every gap it overlaps, by range and then by gap (the order of
# findOverlapPairs). Returns the index of the range each piece comes
# from ('q'), the pieces' starts and widths, and the seqinfo pintersect()
# leaves (the merge with gap's, on y's seqlevels). NULL where the join
# could differ: a class other than GRanges, metadata columns on 'gap', a
# circular sequence, a seqlength the merge fills in, an empty range, gaps
# overlapping or out of order, or too many seqlevels for exact keys.
#
# @param y the input GRanges, whose class and seqinfo the moved ranges have
# @param gap the gaps of the excluded regions, on strand "*"
# @param x_code the index of each range's sequence in seqlevels(y)
# @param x_start,x_end the start and end of each range
intersect_gaps <- function(y, gap, x_code, x_start, x_end) {
  si <- merge(seqinfo(y), seqinfo(gap))
  lev <- seqlevels(y)
  if (!identical(as.vector(class(y)), "GRanges") ||
      !identical(as.vector(class(gap)), "GRanges") ||
      length(mcols(gap)) != 0L || any(Seqinfo::isCircular(si) %in% TRUE) ||
      length(lev) >= 2097152L) {
    return(NULL)
  }
  # the seqinfo pintersect() gives the result, made on an empty copy
  x0 <- y[0L]
  seqinfo(x0) <- si
  Seqinfo::seqlevels(x0) <- lev
  if (!identical(seqlengths(x0), seqlengths(y))) {
    return(NULL)
  }
  gap_code <- match(seqlevels(gap), lev)[as.integer(seqnames(gap))]
  on_x <- which(!is.na(gap_code))
  gap_code <- gap_code[on_x]
  # stable, so each sequence's gaps stay in index order
  o <- order(gap_code, method = "radix")
  gap_code <- gap_code[o]
  gap_start <- start(gap)[on_x][o]
  gap_end <- end(gap)[on_x][o]
  n_gap <- length(gap_code)
  if (any(x_end < x_start) || any(gap_end < gap_start) ||
      any(gap_code[-1L] == gap_code[-n_gap] &
          gap_start[-1L] <= gap_end[-n_gap])) {
    return(NULL)
  }
  # positions offset by 2^32 per sequence, so the gaps of all sequences
  # form one sorted vector and no key crosses into another sequence
  x_off <- x_code * 4294967296
  gap_off <- gap_code * 4294967296
  first <- findInterval(x_off + x_start, gap_off + gap_end,
                        left.open = TRUE) + 1L
  last <- findInterval(x_off + x_end, gap_off + gap_start)
  n_hit <- pmax.int(last - first + 1L, 0L)
  q <- rep.int(seq_along(x_start), n_hit)
  g <- sequence(n_hit, from = first)
  new_start <- pmax.int(x_start[q], gap_start[g])
  new_end <- pmin.int(x_end[q], gap_end[g])
  list(q = q, start = new_start, width = new_end - new_start + 1L,
       seqinfo = seqinfo(x0))
}

# The hits of findOverlaps(bait blocks of replicate r, y), for every r in
# one pass, in the order of the separate calls: by bait block, and within
# one block by index in 'y' where IRanges builds its NCList on the bait
# blocks (when replicate r has fewer of them on that sequence than 'y'
# has ranges), in the order of y's NCList otherwise.
#
# @param bait the bait blocks of all replicates, in replicate order
# @param bait_code the index of each block's sequence in seqlevels(y)
# @param bait_rep the replicate of each block
# @param n_rep the number of replicates
# @param y the input GRanges
bait_overlaps <- function(bait, bait_code, bait_rep, n_rep, y) {
  n_lev <- length(seqlevels(y))
  rep_chr <- (bait_rep - 1L) * n_lev + bait_code
  n_bait <- tabulate(rep_chr, n_rep * n_lev)[rep_chr]
  n_y <- tabulate(as.integer(seqnames(y)), n_lev)[bait_code]
  on_bait <- which(n_bait < n_y)
  on_y <- which(n_bait >= n_y)
  # a GNCList query puts the NCList on the bait blocks; on_y's blocks
  # outnumber y's ranges on each sequence, so IRanges puts it on 'y'
  fo_bait <- if (length(on_y) == 0L) {
    findOverlaps(GNCList(bait), y)
  } else {
    findOverlaps(GNCList(bait[on_bait]), y)
  }
  query <- on_bait[queryHits(fo_bait)]
  subject <- subjectHits(fo_bait)
  if (length(on_y) != 0L) {
    fo_y <- findOverlaps(bait[on_y], y)
    query <- c(query, on_y[queryHits(fo_y)])
    subject <- c(subject, subjectHits(fo_y))
    # stable: each block keeps its hits' order
    o <- order(query, method = "radix")
    query <- query[o]
    subject <- subject[o]
  }
  list(query = query, subject = subject)
}

# Segmentation block bootstrap sub-function: the random bait blocks of
# 'n_rep' replicates, drawn in the order one replicate at a time drew
# them, and the blocks tiling the segmentation they move to
#
# @param n_rep the number of replicates
# @param y the input GRanges
# @param seg the segmentation ranges
# @param state the segmentation state (along seg)
# @param L_b the length of the blocks
# @param chr_names the chromosomes of the segmentation ranges
# @param chr_lens the chromosome lengths
# @param proportion_length whether to use scaled block length
#        (scaling by the proportion of the segmentation state
#         out of the total genome length)
#
# @return a list: 'draws', per replicate the bait block starts ('start')
# and chromosomes ('chr'); 'width', the bait block widths; 'target_start'
# and 'target_chr', the blocks the bait blocks move to
seg_bootstrap_draws <- function(n_rep, y, seg, state, L_b,
                                chr_names, chr_lens,
                                proportion_length = TRUE) {

  # sequence along the states: 1,2,..,S
  esses <- seq_len(max(state))
  chr_names <- as.character(chr_names)

  # break apart this into two sub-functions, defined below
  if (proportion_length) {
    out <- seg_bootstrap_prop(seg, state, L_b, chr_names, chr_lens,
                              esses, n_rep)
  } else {
    out <- seg_bootstrap_noprop(seg, state, L_b, chr_names, chr_lens,
                                esses, n_rep)
  }
  out$target_chr <- factor(out$target_chr, levels = seqlevels(y))
  out
}

# breaking up the seg_bootstrap_draws() code... proportional length case
seg_bootstrap_prop <- function(seg, state, L_b, chr_names, chr_lens,
                               esses, n_rep) {
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

  draws <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
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
      # the tiling blocks do not depend on the draws: computed once, at
      # the point of the first replicate where they were computed before
      if (r > 1L) {
        return(list(random_start_s, random_chr_names_s))
      }
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

    if (r == 1L) {
      # number of blocks for each segmentation state
      nblocks_per_state <- lengths(lapply(obj, `[[`, 1))
      rearranged_blocks_start <- unlist(lapply(obj, `[[`, 3))
      rearr_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 4))
    }

    # unlist the results in 'obj'
    draws[[r]] <- list(
      start = unlist(lapply(obj, `[[`, 1)),
      chr = do.call(c, lapply(obj, `[[`, 2))
    )
  }

  list(
    draws = draws,
    width = rep(L_b_per_state, nblocks_per_state),
    target_start = rearranged_blocks_start,
    target_chr = rearr_blocks_chr_names
  )
}

# breaking up the seg_bootstrap_draws() code... non-proportional length case
seg_bootstrap_noprop <- function(seg, state, L_b, chr_names, chr_lens,
                                 esses, n_rep) {

  # number of rearranged blocks within each range of the segmentation
  nblocks_per_seg <- ceiling(width(seg) / L_b)
  # total number of tiling blocks
  n <- sum(nblocks_per_seg)

  draws <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    # sample segments according to segmentation length
    random_seg <- sample(length(seg), n, replace = TRUE, prob = width(seg))
    # random starts within segments, as specified by 'nblocks_per_seg'
    rmin <- start(seg)[random_seg]
    rmax <- (nblocks_per_seg[random_seg] - 1) * L_b + rmin
    random_start <- round(runif(n, rmin, rmax))
    if (r == 1L) {
      # start positions of rearranged blocks tiling the segments, which
      # do not depend on the draws
      rearr_blocks_start_list <- lapply(seq_along(seg), function(j) {
        seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
      })
      # start positions and chromosomes of blocks tiling segments of
      # state 's'
      rearr_blocks_start <- lapply(esses, function(s) {
        unlist(rearr_blocks_start_list[state == s], use.names = FALSE)
      })
      rearr_blocks_chr_names <- lapply(esses, function(s) {
        rep(chr_names[state == s], nblocks_per_seg[state == s])
      })
    }

    # loop over states 's', splitting out the random segments for 's'
    obj <- lapply(esses, function(s) {
      # index of those random segments chosen that are state 's'
      idx <- which(state[random_seg] == s)
      # random starts for state s
      random_start_s <- random_start[idx]
      # names of chromosomes for random segments that are state 's'
      random_chr_names_s <- chr_names[random_seg[idx]]
      # see if we have too few or too many blocks for this state
      rand_n <- length(random_start_s)
      rearr_n <- length(rearr_blocks_start[[s]])
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
        random_chr_names_s[rand_idx]
      ))
    })

    # unlist the results in 'obj'
    draws[[r]] <- list(
      start = unlist(lapply(obj, `[[`, 1)),
      chr = do.call(c, lapply(obj, `[[`, 2))
    )
  }

  list(
    draws = draws,
    width = L_b,
    target_start = unlist(rearr_blocks_start),
    target_chr = do.call(c, rearr_blocks_chr_names)
  )
}
