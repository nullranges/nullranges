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
#' @useDynLib nullranges, .registration = TRUE, .fixes = "C_"
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
        n_rep, unseg_setup$L_b, unseg_setup$chrom_lens, unseg_setup$blocks,
        compiled = TRUE
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
      br <- lapply(bait_draws(blocks), function(d) {
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

  # from its runs, without a value per range
  keep <- lens > 0L
  mcols(br)$iter <- Rle(factor(seq_len(R), levels=seq_len(R))[keep], lens[keep])
  if (storeBlockLength) {
    mcols(br)$blockLength <- Rle(as.integer(blockLength))
  }
  mcols(br)$block <- NULL # this is not needed except for debugging perhaps
  # new("BootRanges", br) would run validObject over every range again; a
  # plain GRanges here is already valid and within its seqlengths
  if (!identical(as.vector(class(br)), "GRanges")) {
    return(new("BootRanges", br))
  }
  out <- new("BootRanges")
  for (s in slotNames("GRanges")) {
    slot(out, s) <- slot(br, s)
  }
  out
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
# GRanges, whose 'block' column holds the replicate where the R steps
# below count the rows per replicate from it (bootRanges() drops the
# column). Returns the ranges and the number per replicate, or NULL where
# one replicate at a time could differ (a replicate with fewer bait
# blocks than targets, a bait block on a sequence 'y' lacks or past the
# end of a circular one, an input shift_swap_trim() declines). The moved positions are computed on
# vectors, by moved_positions() where it can, with the exclude step when
# exclude_tables() prepares it, and 'y' is subset once, for the ranges kept.
move_blocks_bulk <- function(y, blocks, exclude, excludeOption) {
  n_blocks <- length(blocks$target_start)
  if (is.null(blocks$bait_start)) {
    draws <- blocks$draws
    n_rep <- length(draws)
    bait_start <- lapply(draws, `[[`, "start")
    if (any(lengths(bait_start) != n_blocks)) {
      return(NULL)
    }
    bait_start <- unlist(bait_start, use.names = FALSE)
    bait_chr <- unlist(lapply(draws, `[[`, "chr"), use.names = FALSE)
    bait_code <- match(bait_chr, seqlevels(y))
  } else {
    n_rep <- length(blocks$bait_len)
    if (any(blocks$bait_len != n_blocks)) {
      return(NULL)
    }
    bait_start <- blocks$bait_start
    bait_chr <- NULL
    bait_code <- match(blocks$chr_names, seqlevels(y))[blocks$bait_seg]
  }
  circular <- unname(Seqinfo::isCircular(y)) %in% TRUE
  if (anyNA(bait_code) || any(circular[bait_code] & (bait_start < 1 |
      bait_start > unname(seqlengths(y))[bait_code]))) {
    return(NULL)
  }
  bait_ranges <- IRanges(
    start = bait_start,
    width = rep_len(blocks$width, length(bait_start))
  )
  block_shift <- rep(blocks$target_start, n_rep) - start(bait_ranges)
  target_chr <- blocks$target_chr
  if (is(target_chr, "Rle")) {
    target_chr <- decode(target_chr)
  }
  # a warning or error here is raised again by the steps below, where the
  # exclude step runs without these tables. Without 'exclude', the tables
  # of no exclude step, so that the rows go to moved_ranges() as
  # moved_positions() gives them, with the metadata columns it gathers.
  ex <- NULL
  if (!is.null(exclude)) {
    ex <- tryCatch(
      exclude_tables(y, exclude, excludeOption),
      warning = function(w) NULL, error = function(e) NULL
    )
  } else {
    ex <- list(mode = 0L, seg = integer(), start = integer(),
               end = integer(), seqinfo = NULL)
  }
  moved <- moved_positions(
    y, bait_ranges, bait_code, n_rep, block_shift, target_chr, ex
  )
  if (isFALSE(moved)) {
    ex <- NULL
    moved <- moved_positions(
      y, bait_ranges, bait_code, n_rep, block_shift, target_chr
    )
  }
  if (is.null(moved)) {
    ex <- NULL
    # moved_positions() reads only starts and widths; with no NA in
    # bait_code this GRanges cannot fail, so building it here changes no
    # fallback
    if (is.null(bait_chr)) {
      bait_chr <- blocks$chr_names[blocks$bait_seg]
    }
    random_blocks <- GRanges(seqnames = bait_chr, ranges = bait_ranges)
    bait_rep <- rep(seq_len(n_rep), each = n_blocks)
    fo <- bait_overlaps(random_blocks, bait_code, bait_rep, n_rep, y)
    hit_block <- fo$query
    moved <- shift_swap_trim_positions(
      y, fo$subject, block_shift[hit_block], target_chr, n_rep, hit_block
    )
    if (is.null(moved)) {
      return(NULL)
    }
    i <- which(moved$width > 0L)
    moved <- list(
      s = moved$s[i], code = moved$code[i], start = moved$start[i],
      width = moved$width[i], block = bait_rep[hit_block[i]],
      chr_attr = moved$chr_attr
    )
  }
  # only rows of positive width from here: overlapsAny() is per range
  s <- moved$s
  code <- moved$code
  new_start <- moved$start
  new_width <- moved$width
  # moved_positions() gives the rows per replicate ('lens') in place of
  # each row's replicate ('block'), which the R exclude steps below label
  # their rows with
  block <- moved$block
  lens <- moved$lens
  if (is.null(ex) && !is.null(exclude) && !is.null(lens)) {
    block <- rep.int(seq_len(n_rep), lens)
    lens <- NULL
  }
  if (!is.null(ex)) {
    # the exclude step, if any, is applied: the rows kept, or the pieces
    y_prime <- moved_ranges(
      y, s, code, new_start, new_width, moved$chr_attr, NULL, ex$seqinfo,
      moved$run_code, moved$run_length, moved$cols, with_mcols = TRUE
    )
  } else if (!is.null(exclude) && excludeOption == "drop") {
    i <- outside_exclude(y, exclude, code, new_start, new_width)
    if (is.null(i)) {
      i <- which(!overlapsAny(
        moved_ranges(y, s, code, new_start, new_width, moved$chr_attr),
        exclude, type = "any"
      ))
    }
    y_prime <- moved_ranges(
      y, s[i], code[i], new_start[i], new_width[i], moved$chr_attr, block[i]
    )
  } else if (!is.null(exclude) && excludeOption == "trim") {
    gap <- gaps(exclude, end = seqlengths(y))
    gap <- gap[strand(gap) == "*"]
    clip <- intersect_gaps(
      y, gap, code, new_start, new_start + (new_width - 1L)
    )
    if (is.null(clip)) {
      y_prime <- moved_ranges(
        y, s, code, new_start, new_width, moved$chr_attr, block
      )
      y_prime <- plyranges::join_overlap_intersect(y_prime, gap)
    } else {
      i <- clip$q
      y_prime <- moved_ranges(
        y, s[i], code[i], clip$start, clip$width, moved$chr_attr, block[i],
        clip$seqinfo
      )
    }
  } else {
    y_prime <- moved_ranges(
      y, s, code, new_start, new_width, moved$chr_attr, block, NULL,
      moved$run_code, moved$run_length, with_mcols = TRUE
    )
  }
  if (is.null(lens)) {
    lens <- tabulate(mcols(y_prime)$block, n_rep)
  }
  list(ranges = y_prime, lens = lens)
}

# What shift_swap_trim_positions() gives for the hits of bait_overlaps(),
# in one compiled pass, restricted to the rows of positive width, with
# the number of rows of each replicate as 'lens' (the rows are in
# replicate order): for each bait block in order, the ranges of 'y' it
# overlaps by increasing index in 'y', found as sorted_overlaps() finds
# them, and the runs of the code column
# ('run_code', 'run_length'). NULL wherever bait_overlaps() would not
# take sorted_overlaps() for every block, either function could decline,
# or there is no hit: the caller then runs them. With 'ex', the rows
# outside_exclude() keeps (drop) or the pieces intersect_gaps() cuts (trim),
# in the same order; FALSE where that step declines. Given 'ex', also the
# metadata columns of 'y' at the rows' indices in 'y', as 'cols', in place
# of 's', when gather_columns() gives them.
#
# @param y the input GRanges
# @param bait the ranges (IRanges) of the bait blocks of all replicates,
# in replicate order
# @param bait_code the index of each block's sequence in seqlevels(y)
# @param n_rep the number of replicates
# @param block_shift the shift of each block
# @param target_chr the sequence each block moves to, a factor
# @param ex the output of exclude_tables(), the tables of no exclude step
# (mode 0), or NULL
moved_positions <- function(y, bait, bait_code, n_rep, block_shift,
                            target_chr, ex = NULL) {
  chr_attr <- attributes(target_chr)
  lev <- seqlevels(y)
  if (!identical(as.vector(class(y)), "GRanges") ||
      !is.numeric(block_shift) || length(lev) >= 2097152L ||
      !is.factor(target_chr) ||
      !setequal(names(chr_attr), c("levels", "class")) ||
      !identical(levels(target_chr), lev)) {
    return(NULL)
  }
  circular <- unname(Seqinfo::isCircular(y)) %in% TRUE
  seqlen <- unname(seqlengths(y))
  cols <- NULL
  if (is.null(ex)) {
    ex <- list(mode = 0L, seg = integer(), start = integer(), end = integer())
  } else {
    cols <- gather_columns(y)
  }
  moved <- .Call(
    C_moved_rows_lens, as.integer(bait_code), start(bait), width(bait),
    as.integer(n_rep), as.integer(seqnames(y)), start(y), end(y), width(y),
    circular, as.integer(seqlen), !(circular | is.na(seqlen)),
    suppressWarnings(as.integer(block_shift)), as.integer(target_chr),
    ex$mode, ex$seg, ex$start, ex$end, cols
  )
  if (!is.list(moved)) {
    return(moved)
  }
  moved$chr_attr <- chr_attr
  moved
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
# metadata columns at 's' (none and no names unless 'with_mcols'), with
# 'block' added as a column when given, the given starts, widths and
# seqlevel codes, and 'seqinfo' (y's when NULL). A strand of one run is
# not subset. Given the runs of 'code', maximal, the seqnames
# are built from them.
#
# @param y the input GRanges
# @param s the index in 'y' of each range
# @param code the index of each range's sequence in seqlevels(y)
# @param start,width the start and width of each range
# @param chr_attr the attributes of the seqnames factor
# @param block the replicate of each range, or NULL
# @param seqinfo the seqinfo of the result, or NULL
# @param run_code,run_length the runs of 'code', or NULL
# @param cols the metadata columns of 'y' at 's', from moved_positions(),
# or NULL; given them, 's' is not read
# @param with_mcols whether to keep y's names and metadata columns
moved_ranges <- function(y, s, code, start, width, chr_attr, block = NULL,
                         seqinfo = NULL, run_code = NULL,
                         run_length = NULL, cols = NULL,
                         with_mcols = !is.null(block)) {
  if (is.null(run_length)) {
    attributes(code) <- chr_attr
    ans_seqnames <- Rle(code)
  } else {
    attributes(run_code) <- chr_attr
    ans_seqnames <- Rle(run_code, run_length)
  }
  n <- length(if (is.null(cols)) s else start)
  ans_strand <- strand(y)
  one_strand <- nrun(ans_strand) == 1L
  plain <- one_strand && with_mcols && plain_rows(y, start, width, n)
  if (!is.null(cols) && !plain) {
    stop("internal error: gathered metadata columns off the plain path")
  }
  if (!one_strand || (with_mcols && !plain)) {
    # checked once for every slot subset, as y[s] does
    s <- normalizeSingleBracketSubscript(s, y, as.NSBS = TRUE)
  }
  ans_strand <- if (one_strand) {
    Rle(runValue(ans_strand), n)
  } else {
    extractROWS(ans_strand, s)
  }
  if (!with_mcols) {
    ans_ranges <- new2("IRanges", start = start, width = width, check = FALSE)
    ans_mcols <- new2("DFrame", nrows = n, check = FALSE)
  } else if (plain) {
    # what the branch below gives, without gathering starts and widths
    # that update_ranges() replaces; 's' holds indices in 1..length(y)
    ans_ranges <- BiocGenerics:::replaceSlots(
      ranges(y), start = start, width = width, check = FALSE
    )
    ans_mcols <- mcols(y, use.names = FALSE)
    if (is.null(cols)) {
      cols <- lapply(as.list(ans_mcols), `[`, s)
    }
    ans_mcols <- BiocGenerics:::replaceSlots(
      ans_mcols, listData = cols, nrows = n, check = FALSE
    )
  } else {
    ans_ranges <- update_ranges(
      extractROWS(ranges(y), s), start = start, width = width
    )
    ans_mcols <- extractROWS(mcols(y, use.names = FALSE), s)
  }
  if (!is.null(block)) {
    ans_mcols$block <- block
  }
  ans <- update(y, seqnames = ans_seqnames, ranges = ans_ranges,
                strand = ans_strand, elementMetadata = ans_mcols,
                check = FALSE)
  if (!is.null(seqinfo)) {
    ans <- update(ans, seqinfo = seqinfo, check = FALSE)
  }
  ans
}

# TRUE when y[s] with the given starts and widths is the ranges of 'y'
# with those starts and widths and each metadata column subset by `[`:
# plain IRanges without names or metadata columns, a DFrame without row
# names or column metadata whose columns are atomic vectors without
# attributes, and integer starts and widths of length 'n', which
# update_ranges() would accept.
#
# @param y the input GRanges
# @param start,width the start and width of each range
# @param n the number of ranges
plain_rows <- function(y, start, width, n) {
  is.integer(start) && length(start) == n &&
    is.integer(width) && length(width) == n &&
    plain_ranges_mcols(y)
}

# The ranges and metadata columns part of plain_rows()
#
# @param y the input GRanges
plain_ranges_mcols <- function(y) {
  rng <- y@ranges
  mc <- y@elementMetadata
  identical(as.vector(class(rng)), "IRanges") &&
    is.null(rng@NAMES) && is.null(rng@elementMetadata) &&
    identical(as.vector(class(mc)), "DFrame") &&
    is.null(mc@rownames) && is.null(mc@elementMetadata) &&
    all(vapply(mc@listData, function(v) {
      is.atomic(v) && is.null(attributes(v))
    }, NA))
}

# The metadata columns of 'y' (its DFrame's listData) when moved_ranges()
# takes its plain branch for them and each is a double, integer or logical
# vector of y's length without attributes, which the compiled pass gathers
# bit for bit as `[` does; NULL otherwise.
#
# @param y the input GRanges
gather_columns <- function(y) {
  if (!identical(as.vector(class(y)), "GRanges") ||
      nrun(y@strand) != 1L || !plain_ranges_mcols(y)) {
    return(NULL)
  }
  cols <- y@elementMetadata@listData
  n <- length(y)
  if (!all(vapply(cols, function(v) {
    typeof(v) %in% c("double", "integer", "logical") && length(v) == n
  }, NA))) {
    return(NULL)
  }
  cols
}

# which(!overlapsAny(x, exclude)) for the moved ranges 'x' of positive
# width whose seqlevel codes, starts and widths are given, on y's seqinfo:
# 'exclude' merged per sequence into sorted disjoint intervals, then, in
# C, a range overlaps iff the first interval ending at or after its start
# begins at or before its end. merge() of the seqinfo runs first, for the
# warning or error findOverlaps() gives. NULL where overlapsAny() could
# differ: an 'exclude' that is not a plain GRanges, a strand other than
# "*", a zero-width range, a circular sequence, or too many seqlevels for
# exact keys.
#
# @param y the input GRanges, whose seqinfo the moved ranges have
# @param exclude the excluded regions
# @param code the index of each range's sequence in seqlevels(y)
# @param start,width the start and width of each range, width > 0
outside_exclude <- function(y, exclude, code, start, width) {
  if (!is.integer(code) || !is.integer(start) || !is.integer(width)) {
    return(NULL)
  }
  iv <- exclude_intervals(y, exclude)
  if (is.null(iv)) {
    return(NULL)
  }
  .Call(C_outside_intervals, code, start, width, iv$seg, iv$start, iv$end)
}

# The intervals outside_exclude() searches: 'exclude' merged per seqlevel
# of 'y' into sorted disjoint intervals, sequence c holding seg[c] + 1 to
# seg[c + 1] of 'start' and 'end'. NULL where outside_exclude() declines
# for any moved ranges.
#
# @param y the input GRanges
# @param exclude the excluded regions
exclude_intervals <- function(y, exclude) {
  if (!identical(as.vector(class(exclude)), "GRanges")) {
    return(NULL)
  }
  si <- merge(seqinfo(y), seqinfo(exclude))
  lev <- seqlevels(y)
  if (length(lev) >= 2097152L || any(Seqinfo::isCircular(si) %in% TRUE) ||
      !all(strand(exclude) == "*") || any(width(exclude) == 0L)) {
    return(NULL)
  }
  ex_code <- match(seqlevels(exclude), lev)[as.integer(seqnames(exclude))]
  on_y <- which(!is.na(ex_code))
  ex_code <- ex_code[on_y]
  n_ex <- length(on_y)
  if (n_ex == 0L) {
    int_start <- int_end <- integer()
    seg <- integer(length(lev) + 1L)
  } else {
    # positions offset by 2^32 per sequence, so the intervals of all
    # sequences form one sorted vector and cummax() stays in a sequence
    key <- ex_code * 4294967296
    ex_start <- key + start(exclude)[on_y]
    o <- order(ex_start, method = "radix")
    ex_start <- ex_start[o]
    ex_end <- cummax((key[o] + end(exclude)[on_y][o]))
    first <- c(TRUE, ex_start[-1L] > ex_end[-n_ex])
    last <- c(which(first)[-1L] - 1L, n_ex)
    int_code <- ex_code[o][first]
    int_start <- as.integer(ex_start[first] - int_code * 4294967296)
    int_end <- as.integer(ex_end[last] - int_code * 4294967296)
    seg <- c(0L, cumsum(tabulate(int_code, length(lev))))
  }
  list(seg = seg, start = int_start, end = int_end)
}

# What moved_positions() needs to apply the exclude step itself: for
# "drop", exclude_intervals(); for "trim", the gaps of 'exclude' as
# gap_tables() orders them, and the seqinfo of the result. The mode is 1
# for drop and 2 for trim. NULL where outside_exclude() or intersect_gaps()
# declines for any moved ranges of positive width, or for another option.
#
# @param y the input GRanges
# @param exclude the excluded regions
# @param excludeOption "drop" or "trim"
exclude_tables <- function(y, exclude, excludeOption) {
  if (excludeOption == "drop") {
    iv <- exclude_intervals(y, exclude)
    if (is.null(iv)) {
      return(NULL)
    }
    return(list(mode = 1L, seg = iv$seg, start = iv$start, end = iv$end,
                seqinfo = NULL))
  }
  if (excludeOption == "trim") {
    gap <- gaps(exclude, end = seqlengths(y))
    gap <- gap[strand(gap) == "*"]
    tb <- gap_tables(y, gap)
    if (is.null(tb) || !tb$ordered) {
      return(NULL)
    }
    return(list(mode = 2L, seg = tb$seg, start = tb$start, end = tb$end,
                seqinfo = tb$seqinfo))
  }
  NULL
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
  tb <- gap_tables(y, gap)
  if (is.null(tb) || any(x_end < x_start) || !tb$ordered) {
    return(NULL)
  }
  gap_code <- tb$code
  gap_start <- tb$start
  gap_end <- tb$end
  clip <- NULL
  if (is.null(attributes(x_code)) && is.null(attributes(x_start)) &&
      is.null(attributes(x_end))) {
    clip <- .Call(
      C_intersect_gaps, x_code, x_start, x_end, tb$seg, gap_start, gap_end
    )
  }
  if (is.null(clip)) {
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
    clip <- list(q, new_start, new_end - new_start + 1L)
  }
  list(q = clip[[1L]], start = clip[[2L]], width = clip[[3L]],
       seqinfo = tb$seqinfo)
}

# The gaps intersect_gaps() clips to, on y's seqlevels: their codes,
# starts and ends ordered by sequence and then index, 'seg' (sequence c
# holds gaps seg[c] + 1 to seg[c + 1]), whether they are disjoint and in
# order along each sequence ('ordered'), and the seqinfo of the result.
# NULL where intersect_gaps() declines whatever the ranges.
#
# @param y the input GRanges, whose class and seqinfo the moved ranges have
# @param gap the gaps of the excluded regions, on strand "*"
gap_tables <- function(y, gap) {
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
  list(
    code = gap_code, start = gap_start, end = gap_end,
    seg = c(0L, cumsum(tabulate(gap_code, length(lev)))),
    ordered = !(any(gap_end < gap_start) ||
      any(gap_code[-1L] == gap_code[-n_gap] &
          gap_start[-1L] <= gap_end[-n_gap])),
    seqinfo = seqinfo(x0)
  )
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
  fo_bait <- sorted_overlaps(bait, bait_code, on_bait, y)
  if (is.null(fo_bait)) {
    # a GNCList query puts the NCList on the bait blocks; on_y's blocks
    # outnumber y's ranges on each sequence, so IRanges puts it on 'y'
    fo_bait <- if (length(on_y) == 0L) {
      findOverlaps(GNCList(bait), y)
    } else {
      findOverlaps(GNCList(bait[on_bait]), y)
    }
    fo_bait <- list(
      query = on_bait[queryHits(fo_bait)], subject = subjectHits(fo_bait)
    )
  }
  query <- fo_bait$query
  subject <- fo_bait$subject
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

# What findOverlaps(GNCList(bait[on_bait]), y) returns, as the blocks'
# indices in 'bait' ('query') and the hits' in 'y' ('subject'): by block,
# and within a block by increasing index in 'y'. Within each sequence the
# ranges of 'y' with start <= the block's end form a prefix in index
# order, and those whose running maximum of ends is below the block's
# start form a prefix of that; between the two, only a range ending
# before an earlier one's end can miss the block. NULL where the call
# could differ: no blocks, a class other than GRanges, a circular
# sequence under a block, an empty range, starts of 'y' out of order
# within a sequence, or too many seqlevels for exact keys.
#
# @param bait the bait blocks of all replicates
# @param bait_code the index of each block's sequence in seqlevels(y)
# @param on_bait the blocks to find the hits of
# @param y the input GRanges
sorted_overlaps <- function(bait, bait_code, on_bait, y) {
  if (length(on_bait) == 0L || !identical(as.vector(class(y)), "GRanges") ||
      length(seqlevels(y)) >= 2097152L) {
    return(NULL)
  }
  b_code <- bait_code[on_bait]
  if (any(unname(Seqinfo::isCircular(y))[b_code] %in% TRUE)) {
    return(NULL)
  }
  b_start <- start(bait)[on_bait]
  b_end <- end(bait)[on_bait]
  y_code <- as.integer(seqnames(y))
  y_start <- start(y)
  y_end <- end(y)
  if (any(b_end < b_start) || any(y_end < y_start)) {
    return(NULL)
  }
  o <- NULL
  if (is.unsorted(y_code)) {
    o <- order(y_code, method = "radix")
    y_code <- y_code[o]
    y_start <- y_start[o]
    y_end <- y_end[o]
  }
  # positions offset by 2^32 per sequence: one vector, sorted when each
  # sequence's part is
  y_off <- y_code * 4294967296
  if (is.unsorted(y_off + y_start)) {
    return(NULL)
  }
  end_max <- cummax(y_off + y_end) - y_off
  last <- count_at_or_before(y_code, y_start, b_code, b_end)
  first <- count_at_or_before(y_code, end_max, b_code, b_start - 1) + 1L
  n_hit <- pmax.int(last - first + 1L, 0L)
  subject <- sequence(n_hit, from = first)
  short <- y_end < end_max
  if (any(short)) {
    n_short <- c(0L, cumsum(short))
    flag <- which(n_hit > 0L & n_short[last + 1L] > n_short[first])
    if (length(flag) != 0L) {
      at <- sequence(n_hit[flag], from = (cumsum(n_hit) - n_hit)[flag] + 1L)
      miss <- y_end[subject[at]] < rep.int(b_start[flag], n_hit[flag])
      if (any(miss)) {
        subject <- subject[-at[miss]]
        n_hit <- n_hit - tabulate(
          rep.int(flag, n_hit[flag])[miss], length(n_hit)
        )
      }
    }
  }
  if (!is.null(o)) {
    subject <- o[subject]
  }
  list(query = rep.int(on_bait, n_hit), subject = subject)
}

# For each query, the number of entries of a list sorted by sequence and
# then position that come before it or at its position on its sequence:
# findInterval() on keys code * 2^32 + position, from a table of
# fixed-width buckets along each sequence (with an empty one before and
# after each sequence's entries) and a step per entry in the query's
# bucket; findInterval() itself for queries in crowded buckets.
#
# @param code,pos the entries' sequences and positions, sorted
# @param x_code,x_pos the queries' sequences, each holding an entry, and
#        positions
count_at_or_before <- function(code, pos, x_code, x_pos) {
  n <- length(code)
  lev <- seq_len(code[n])
  lo <- as.numeric(pos[match(lev, code)])
  hi <- as.numeric(pos[n + 1L - match(lev, rev(code))])
  width <- 2^max(0, ceiling(log2(sum(hi - lo + 1, na.rm = TRUE) / (8 * n))))
  n_bucket <- floor((hi - lo) / width) + 3
  n_bucket[is.na(n_bucket)] <- 0
  base <- cumsum(n_bucket) - n_bucket
  counts <- tabulate(
    base[code] + floor((pos - lo[code]) / width) + 2, sum(n_bucket)
  )
  to <- cumsum(counts)
  from <- to - counts
  j <- base[x_code] +
    pmax(pmin(floor((x_pos - lo[x_code]) / width) + 2, n_bucket[x_code]), 1)
  cnt <- from[j]
  to <- to[j]
  i <- which(to > cnt)
  for (step in 1:16) {
    if (length(i) == 0L) {
      break
    }
    i <- i[pos[cnt[i] + 1L] <= x_pos[i]]
    cnt[i] <- cnt[i] + 1L
    i <- i[cnt[i] < to[i]]
  }
  if (length(i) != 0L) {
    cnt[i] <- findInterval(
      x_code[i] * 4294967296 + x_pos[i], code * 4294967296 + pos
    )
  }
  cnt
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
# and 'target_chr', the blocks the bait blocks move to. When
# seg_draws_compiled() made the draws, 'draws' is NULL and 'bait_start',
# 'bait_seg' (the index in seg) and 'bait_len' (the number per replicate)
# hold them, with 'chr_names'; bait_draws() gives the list
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
  if (!is.null(out$bait_start)) {
    out$chr_names <- chr_names
  }
  out
}

# blocks$draws, built from the vectors of seg_draws_compiled() if need be
bait_draws <- function(blocks) {
  if (is.null(blocks$bait_start)) {
    return(blocks$draws)
  }
  len <- blocks$bait_len
  first <- cumsum(len) - len
  lapply(seq_along(len), function(r) {
    i <- first[r] + seq_len(len[r])
    list(start = as.double(blocks$bait_start[i]),
         chr = blocks$chr_names[blocks$bait_seg[i]])
  })
}

# The draws of seg_bootstrap_prop() (no 'st') or seg_bootstrap_noprop()
# for all replicates by one call of C_seg_draws, which makes every draw as
# sample() and round(runif()) make it, in the same order. Each table holds
# segments 'index' with their widths as probabilities, 'n' draws per
# replicate, and runif() bounds rmin = start, rmax = (nblocks - 1) * L + start
# per segment, the values the R loop computes per draw. NULL, before any
# draw, wherever a draw could fail, warn or give a start past the integer
# range, so the R loop runs and fails or warns as before.
#
# @param tabs per table, list(index, width, start, nblocks, L, n)
# @param n_rep the number of replicates
# @param st the state of each segment (0 for none), noprop only
# @param rearr the number of tiling blocks per state, noprop only
seg_draws_compiled <- function(tabs, n_rep, st = NULL, rearr = NULL) {
  if (length(tabs) < 1L || !is.integer(n_rep) || length(n_rep) != 1L ||
      is.na(n_rep) || n_rep < 1L) {
    return(NULL)
  }
  n_draw <- integer(length(tabs))
  rmax <- vector("list", length(tabs))
  for (k in seq_along(tabs)) {
    tb <- tabs[[k]]
    n <- tb$n
    L <- tb$L
    if (length(tb$width) < 1L || !is.integer(tb$width) ||
        !is.integer(tb$start) || !is.numeric(tb$nblocks) ||
        !is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
        n > .Machine$integer.max || n != trunc(n) ||
        !is.numeric(L) || length(L) != 1L ||
        anyNA(tb$width) || any(tb$width < 0L) || !any(tb$width > 0L)) {
      return(NULL)
    }
    hi <- (tb$nblocks - 1) * L + tb$start
    if (length(hi) != length(tb$start) || !all(is.finite(hi)) ||
        any(hi < tb$start) || any(tb$start < -2147483645L) ||
        any(hi > 2147483645)) {
      return(NULL)
    }
    n_draw[k] <- as.integer(n)
    rmax[[k]] <- as.double(hi)
  }
  m <- lengths(lapply(tabs, `[[`, "width"))
  .Call(
    C_seg_draws, n_rep,
    as.double(unlist(lapply(tabs, `[[`, "width"), use.names = FALSE)),
    as.double(unlist(lapply(tabs, `[[`, "start"), use.names = FALSE)),
    unlist(rmax, use.names = FALSE),
    as.integer(unlist(lapply(tabs, `[[`, "index"), use.names = FALSE)),
    c(0L, cumsum(m)), n_draw, st, rearr
  )
}

# breaking up the seg_bootstrap_draws() code... proportional length case
seg_bootstrap_prop <- function(seg, state, L_b, chr_names, chr_lens,
                               esses, n_rep) {
  seg_width <- width(seg)
  seg_start <- start(seg)
  # length of genome
  L_g <- sum(chr_lens)
  # segmentation state total length
  L_s <- unlist(lapply(esses, function(s) {
    sum(seg_width[state == s])
  }), use.names = FALSE)
  # the block width for each segmentation state is scaled
  # down to the segmentation state size, e.g. if segmentation state
  # is half of the genome, then the block width is half of L_b
  L_b_per_state <- round(L_b * L_s / L_g)

  # per state, the segmentation ranges, which no draw changes
  per_state <- lapply(esses, function(s) {
    in_s <- state == s
    width_s <- seg_width[in_s]
    # number of blocks within each range of the segmentation
    nblocks_per_rng <- ceiling(width_s / L_b_per_state[s])
    list(
      width = width_s, start = seg_start[in_s],
      chr = chr_names[in_s], nblocks = nblocks_per_rng,
      # total number of tiling blocks
      n = sum(nblocks_per_rng), index = which(in_s)
    )
  })

  out <- seg_prop_compiled(per_state, L_b, L_b_per_state, n_rep)
  if (!is.null(out)) {
    return(out)
  }

  draws <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    # loops over states 's'
    obj <- lapply(esses, function(s) {
      # the length of the blocks for this state
      L_bs <- L_b_per_state[s]
      p <- per_state[[s]]
      nblocks_per_rng <- p$nblocks
      n <- p$n
      # sample ranges according to their length
      random_rng <- sample(length(p$width), n,
        replace = TRUE,
        prob = p$width
      )
      # random starts within the ranges, as spec by 'nblocks per rng'
      rmin <- p$start[random_rng]
      rmax <- (nblocks_per_rng[random_rng] - 1) * L_bs + rmin
      random_start_s <- round(runif(n, rmin, rmax))
      # names of chromosomes for random ranges of state 's'
      random_chr_names_s <- p$chr[random_rng]
      # the tiling blocks do not depend on the draws: computed once, at
      # the point of the first replicate where they were computed before
      if (r > 1L) {
        return(list(random_start_s, random_chr_names_s))
      }
      # start positions of blocks tiling ranges of state 's'
      end_s <- p$start + p$width - 1L
      rearr_blocks_start_s <- unlist(lapply(seq_along(p$width), function(j) {
        seq(from = p$start[j], to = end_s[j], by = L_bs)
      }), use.names = FALSE)
      # names of chromosomes for blocks tiling ranges of state 's'
      rearr_blocks_chr_names_s <- rep(p$chr, nblocks_per_rng)
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

# unlist(lapply(seq_along(from), function(j) seq(from[j], to[j], by = by)),
# use.names = FALSE) in one pass. With integer from and to and a double by,
# seq.default computes from + (0:n) * by, n = as.integer((to - from)/by +
# 1e-10), capped at to, and returns the integer from when to == from; the
# seq() loop runs wherever it could fail or take another branch, and where
# every segment is one base wide (its value is then integer).
tile_starts <- function(from, to, by) {
  if (is.integer(from) && is.integer(to) && !is.object(from) &&
      !is.object(to) && length(from) > 0L && length(to) == length(from) &&
      is.double(by) && length(by) == 1L && is.null(attributes(by)) &&
      is.finite(by) && by > 0 && !anyNA(from) && !anyNA(to)) {
    del <- as.double(to) - as.double(from)
    n <- del / by
    if (all(del >= 0) && any(del > 0) && all(n <= .Machine$integer.max)) {
      k <- as.integer(n + 1e-10)
      if (sum(as.double(k)) + length(k) <= .Machine$integer.max) {
        cnt <- k + 1L
        return(pmin(rep.int(from, cnt) + (sequence(cnt) - 1L) * by,
                    rep.int(to, cnt)))
      }
    }
  }
  unlist(lapply(seq_along(from), function(j) {
    seq(from = from[j], to = to[j], by = by)
  }), use.names = FALSE)
}

# seg_bootstrap_prop()'s value by seg_draws_compiled(), the tiling blocks
# computed first; NULL where the R loop must run
seg_prop_compiled <- function(per_state, L_b, L_b_per_state, n_rep) {
  if (length(L_b) != 1L) {
    return(NULL)
  }
  tiles <- tryCatch(lapply(seq_along(per_state), function(s) {
    p <- per_state[[s]]
    end_s <- p$start + p$width - 1L
    list(
      tile_starts(p$start, end_s, L_b_per_state[s]),
      rep(p$chr, p$nblocks)
    )
  }), warning = function(w) NULL, error = function(e) NULL)
  if (is.null(tiles)) {
    return(NULL)
  }
  tabs <- lapply(seq_along(per_state), function(s) {
    p <- per_state[[s]]
    list(index = p$index, width = p$width, start = p$start,
         nblocks = p$nblocks, L = L_b_per_state[s], n = p$n)
  })
  d <- seg_draws_compiled(tabs, n_rep)
  if (is.null(d)) {
    return(NULL)
  }
  nblocks_per_state <- vapply(tabs, function(tb) as.integer(tb$n), 1L)
  list(
    draws = NULL,
    width = rep(L_b_per_state, nblocks_per_state),
    target_start = unlist(lapply(tiles, `[[`, 1)),
    target_chr = do.call(c, lapply(tiles, `[[`, 2)),
    bait_start = d[[1L]], bait_seg = d[[2L]], bait_len = d[[3L]]
  )
}

# breaking up the seg_bootstrap_draws() code... non-proportional length case
seg_bootstrap_noprop <- function(seg, state, L_b, chr_names, chr_lens,
                                 esses, n_rep) {

  seg_width <- width(seg)
  seg_start <- start(seg)
  # number of rearranged blocks within each range of the segmentation
  nblocks_per_seg <- ceiling(seg_width / L_b)
  # total number of tiling blocks
  n <- sum(nblocks_per_seg)

  out <- seg_noprop_compiled(seg_width, seg_start, nblocks_per_seg, n,
                             L_b, state, chr_names, esses, n_rep)
  if (!is.null(out)) {
    return(out)
  }

  draws <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    # sample segments according to segmentation length
    random_seg <- sample(length(seg), n, replace = TRUE, prob = seg_width)
    # random starts within segments, as specified by 'nblocks_per_seg'
    rmin <- seg_start[random_seg]
    rmax <- (nblocks_per_seg[random_seg] - 1) * L_b + rmin
    random_start <- round(runif(n, rmin, rmax))
    if (r == 1L) {
      # start positions of rearranged blocks tiling the segments, which
      # do not depend on the draws
      seg_end <- seg_start + seg_width - 1L
      rearr_blocks_start_list <- lapply(seq_along(seg), function(j) {
        seq(from = seg_start[j], to = seg_end[j], by = L_b)
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

    random_state <- state[random_seg]
    # loop over states 's', splitting out the random segments for 's'
    obj <- lapply(esses, function(s) {
      # index of those random segments chosen that are state 's'
      idx <- which(random_state == s)
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

# seg_bootstrap_noprop()'s value by seg_draws_compiled(), the tiling blocks
# computed first; NULL where the R loop must run
seg_noprop_compiled <- function(seg_width, seg_start, nblocks_per_seg, n,
                                L_b, state, chr_names, esses, n_rep) {
  if (length(L_b) != 1L || !is.numeric(state) ||
      length(state) != length(seg_width) || length(esses) < 1L) {
    return(NULL)
  }
  tiles <- tryCatch({
    seg_end <- seg_start + seg_width - 1L
    in_state <- lapply(esses, function(s) which(state == s))
    # segments in no state are tiled too, for seq()'s errors
    none <- rep.int(TRUE, length(seg_width))
    none[unlist(in_state)] <- FALSE
    if (any(none)) {
      tile_starts(seg_start[none], seg_end[none], L_b)
    }
    list(
      lapply(in_state, function(i) {
        tile_starts(seg_start[i], seg_end[i], L_b)
      }),
      lapply(esses, function(s) {
        rep(chr_names[state == s], nblocks_per_seg[state == s])
      })
    )
  }, warning = function(w) NULL, error = function(e) NULL)
  if (is.null(tiles)) {
    return(NULL)
  }
  tabs <- list(list(index = seq_along(seg_width), width = seg_width,
                    start = seg_start, nblocks = nblocks_per_seg, L = L_b,
                    n = n))
  d <- seg_draws_compiled(tabs, n_rep,
                          st = match(state, esses, nomatch = 0L),
                          rearr = lengths(tiles[[1L]]))
  if (is.null(d)) {
    return(NULL)
  }
  list(
    draws = NULL,
    width = L_b,
    target_start = unlist(tiles[[1L]]),
    target_chr = do.call(c, tiles[[2L]]),
    bait_start = d[[1L]], bait_seg = d[[2L]], bait_len = d[[3L]]
  )
}
