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
    # the runs of a GRanges' strand Rle hold every strand present; the Rle
    # comparison makes R's JIT compile S4Vectors' rep method on a second
    # call; the stopifnot text is kept for the failing path
    ex_strand <- if (identical(as.vector(class(exclude)), "GRanges")) strand(exclude)
    if (!(is(ex_strand, "Rle") && isTRUE(all(runValue(ex_strand) == "*")))) {
      stopifnot(all(strand(exclude) == "*"))
    }
  }
  if (!is.null(seg)) {
    seg_mcols <- mcols(seg)
    # the check, with its message, where it fails
    if (!("state" %in% names(seg_mcols))) {
      stopifnot("state" %in% names(mcols(seg)))
    }
  }

  # base's code (bootRanges_base()) where the paths below could differ from
  # it: no blockLength (missing there too, for the same error), R not one
  # number >= 1, a 'y' whose moved ranges base's replicates would reject, a
  # generator .Random.seed may not hold
  if (missing(blockLength)) {
    return(bootRanges_base(y, R = R, seg = seg, exclude = exclude,
                           excludeOption = excludeOption,
                           proportionLength = proportionLength, type = type,
                           withinChrom = withinChrom,
                           storeBlockLength = storeBlockLength))
  }
  if (!plain_count(R) || !boot_y_valid(y) || user_rng()) {
    return(bootRanges_base(y, blockLength, R, seg, exclude, excludeOption,
                           proportionLength, type, withinChrom,
                           storeBlockLength))
  }

  # the arguments base's first replicate reads before its first draw, read
  # in its order and outside the handlers below, which would leave them
  # half evaluated
  if (is.null(seg)) {
    force(type)
    force(blockLength)
    force(withinChrom)
  } else {
    force(proportionLength)
    force(blockLength)
  }
  # a block length other than one whole number: base's code. With a
  # fractional one, its tiling of a segmentation without proportional
  # lengths gives fewer target starts than target sequences, which base
  # pairs per replicate and the moves for all replicates at once do not.
  if (!plain_count(blockLength) || blockLength != floor(blockLength)) {
    return(bootRanges_base(y, blockLength, R, seg, exclude, excludeOption,
                           proportionLength, type, withinChrom,
                           storeBlockLength))
  }
  if (is.null(seg)) {
    # done once for all replicates. Any warning or error it signals, base's
    # first replicate signals before its first draw, or after its draws for
    # the tiling of the bootstrap across chromosomes: base's code runs.
    unseg_setup <- tryCatch(
      unseg_bootstrap_setup(
        y = y, L_b = blockLength, type = type,
        within_chrom = withinChrom
      ),
      warning = function(w) NULL, error = function(e) NULL
    )
    if (is.null(unseg_setup)) {
      return(bootRanges_base(y, blockLength, R, seg, exclude, excludeOption,
                             proportionLength, type, withinChrom,
                             storeBlockLength))
    }
  }

  # as many replicates as replicate(R, ...) makes
  n_rep <- length(integer(R))
  bait <- !is.null(seg) ||
    (unseg_setup$type == "bootstrap" && isFALSE(unseg_setup$within_chrom))

  if (bait) {
    # read before the draws; base reads it after its first replicate is
    # drawn and moved. Where reading it fails, that replicate runs, as base
    # runs it, and then the error is raised.
    ex_error <- tryCatch({
      force(exclude)
      NULL
    }, error = function(e) e)
    if (!is.null(ex_error)) {
      bootRanges_base(y, blockLength, R = 1, seg = seg, exclude = NULL,
                      excludeOption = excludeOption,
                      proportionLength = proportionLength, type = type,
                      withinChrom = withinChrom, storeBlockLength = FALSE)
      stop(ex_error)
    }
    # what moving the blocks reads of 'y' and 'exclude', which no draw
    # changes, so that the compiled draws can move each replicate's blocks
    # while the next are drawn; NULL where they are moved after the draws
    fuse <- tryCatch(
      bulk_setup(y, exclude, excludeOption),
      warning = function(w) NULL, error = function(e) NULL
    )
    # the replicates in chunks of at most boot_block_budget() bait blocks
    # (at least one replicate each), so that the blocks held at once stay
    # bounded; moving draws nothing, so the draws keep their order
    per_rep <- bait_blocks_per_rep(
      seg, seg_mcols, chr_lens, blockLength, proportionLength, unseg_setup
    )
    chunks <- boot_chunks(n_rep, per_rep)
    # where each replicate has at least as many blocks as 'y' has ranges,
    # some sequence holds at least as many of a replicate's blocks as ranges
    # of 'y', and the compiled moves decline it unless y's ranges there are
    # in the order of its NCList (or none): unless they are on every
    # sequence, the draws are made without them and the blocks moved by
    # move_blocks_pieces()
    r_route <- identical(as.vector(class(y)), "GRanges") &&
      !is.na(per_rep) && per_rep >= length(y) && !y_nclist_order(y)
    # each chunk's replicates are drawn before any of them is moved. On any
    # warning or error from here, .Random.seed is put back as it was before
    # the first draw and base's code, which draws and moves one replicate
    # at a time, makes the call, so that the condition comes at its point
    # of the stream.
    mark <- rng_mark()
    moved <- tryCatch({
      # what seg$state gives, without DataFrame's [[ and its selectMethod()
      state <- if (!is.null(seg)) {
        if (identical(as.vector(class(seg_mcols)), "DFrame")) {
          seg_mcols@listData[["state"]]
        } else {
          seg$state
        }
      }
      fz <- if (!r_route) fuse
      # the random bait blocks of n replicates first, in replicate order,
      # then the features they capture moved for those replicates at once
      draw_move <- function(n) {
        blocks <- if (!is.null(seg)) {
          seg_bootstrap_draws(
            n, y = y, seg = ranges(seg),
            state = state, L_b = blockLength,
            chr_names = seqnames(seg),
            chr_lens = chr_lens,
            proportion_length = proportionLength,
            fuse = fz
          )
        } else {
          unseg_across_draws(
            n, unseg_setup$L_b, unseg_setup$chrom_lens, unseg_setup$blocks,
            compiled = TRUE, fuse = fz
          )
        }
        # NULL for any input where the replicates, run one at a time, could
        # warn, fail or differ: then they are run one at a time
        bulk <- tryCatch({
          pieces <- if (r_route) {
            move_blocks_pieces(y, blocks, exclude, excludeOption,
                               boot_block_budget())
          } else {
            FALSE
          }
          if (isFALSE(pieces)) {
            move_blocks_bulk(y, blocks, exclude, excludeOption, fz)
          } else {
            pieces
          }
        }, warning = function(w) NULL, error = function(e) NULL)
        if (!is.null(bulk)) {
          bulk
        } else {
          br <- lapply(bait_draws(blocks), function(d) {
            exclude_ranges(move_blocks(y, blocks, d), y, exclude,
                           excludeOption)
          })
          list(ranges = do.call(c, br), lens = lengths(br))
        }
      }
      if (length(chunks) == 1L) {
        draw_move(n_rep)
      } else {
        bind_chunks(lapply(chunks, draw_move))
      }
    }, warning = function(w) NULL, error = function(e) NULL)
    if (is.null(moved)) {
      rng_restore(mark)
      return(bootRanges_base(y, blockLength, R, seg, exclude, excludeOption,
                             proportionLength, type, withinChrom,
                             storeBlockLength))
    }
    lens <- moved$lens
    br <- moved$ranges
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
  iter <- Rle(factor(seq_len(R), levels=seq_len(R))[keep], lens[keep])
  # storeBlockLength is read once, after the replicates, as base reads it
  if (plain_iter_cols(br, sum(lens)) && isFALSE(storeBlockLength)) {
    # the two mcols<- lines below on the DFrame's columns: list [[<-
    # appends a new name at the end and replaces an existing one in place,
    # as DataFrame's [[<- does
    br_mcols <- br@elementMetadata
    cols <- br_mcols@listData
    cols[["iter"]] <- iter
    cols[["block"]] <- NULL
    br <- BiocGenerics:::replaceSlots(
      br, elementMetadata = BiocGenerics:::replaceSlots(
        br_mcols, listData = cols, check = FALSE
      ),
      check = FALSE
    )
  } else {
    mcols(br)$iter <- iter
    if (storeBlockLength) {
      mcols(br)$blockLength <- Rle(as.integer(blockLength))
    }
    mcols(br)$block <- NULL # this is not needed except for debugging perhaps
  }
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

# TRUE for one finite number >= 1 without a class (for R, integer(R) then
# reads it as replicate(R, ...) in base's code does)
#
# @param x the number of replicates or the block length
plain_count <- function(x) {
  is.numeric(x) && !is.object(x) && length(x) == 1L && is.finite(x) &&
    x >= 1
}

# FALSE where 'y', exactly a GRanges, fails a check of the validity methods
# of GRanges and its superclasses (GenomicRanges 1.65.4, S4Vectors 0.51.10)
# other than the warning for ranges past a seqlength. Base's replicates run
# them on the ranges they move (shift(), seqnames<-); the moves here, on a
# GRanges, do not. TRUE for other classes, whose moves call shift().
#
# @param y the input GRanges
boot_y_valid <- function(y) {
  if (!identical(as.vector(class(y)), "GRanges")) {
    return(TRUE)
  }
  sn <- y@seqnames
  st <- y@strand
  rng <- y@ranges
  mc <- y@elementMetadata
  n <- length(sn)
  nm <- names(rng)
  # GenomicRanges: seqnames, ranges, strand, mcols and seqinfo
  is(sn, "Rle") && is.factor(runValue(sn)) && is.null(names(sn)) &&
    !anyNA(runValue(sn)) &&
    as.vector(class(rng)) %in% c("IRanges", "StitchedIPos", "UnstitchedIPos") &&
    is.factor(runValue(st)) &&
    identical(levels(runValue(st)), c("+", "-", "*")) &&
    !anyNA(runValue(st)) &&
    identical(seqlevels(y@seqinfo), levels(runValue(sn))) &&
    # GRanges: no metadata columns on the ranges, no row names on mcols
    is.null(rng@elementMetadata) &&
    # Vector: slots parallel to 'y', names, mcols a DataFrame
    NROW(rng) == n && NROW(st) == n &&
    (is.null(nm) || (is.character(nm) && is.null(attributes(nm)) &&
      length(nm) == n)) &&
    (is.null(mc) || (is(mc, "DataFrame") && nrow(mc) == n &&
      is.null(rownames(mc)) &&
      !any(c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
             "isCircular", "start", "end", "width", "element") %in%
           colnames(mc))))
}

# TRUE when the generator is user-supplied (RNGkind()), whose state need
# not be in .Random.seed, so that rng_restore() could not rewind it
user_rng <- function() {
  seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  is.integer(seed) && length(seed) >= 1L && !is.na(seed[1L]) &&
    seed[1L] %% 100L == 5L
}

# .Random.seed as it stands, NULL when there is none, for rng_restore()
rng_mark <- function() {
  list(seed = get0(".Random.seed", envir = globalenv(), inherits = FALSE))
}

# Puts back the .Random.seed of rng_mark(): the next draw continues from
# it, or, when there was none, seeds the generator afresh, as the first
# draw of a session does
#
# @param mark the value of rng_mark()
rng_restore <- function(mark) {
  if (!is.null(mark$seed)) {
    assign(".Random.seed", mark$seed, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
}

# TRUE when bootRanges' result 'br' is exactly a GRanges without names
# whose metadata columns are a DFrame without row names, column metadata
# or metadata, with 'n' rows, one per range
#
# @param br the concatenated replicates
# @param n the length of the iter column
plain_iter_cols <- function(br, n) {
  if (!identical(as.vector(class(br)), "GRanges") || !is.null(names(br))) {
    return(FALSE)
  }
  mc <- br@elementMetadata
  identical(as.vector(class(mc)), "DFrame") &&
    is.null(mc@rownames) && is.null(mc@elementMetadata) &&
    length(mc@metadata) == 0L &&
    identical(mc@nrows, length(br)) && isTRUE(n == mc@nrows)
}

# TRUE where, on every sequence, the ranges of 'y' (a GRanges) in index
# order are non-empty and in the order of the NCList IRanges builds on them
# (start ascending, then end descending): the compiled moves then give
# each block's hits in findOverlaps()' order whichever side holds the
# NCList, as moved_setup() checks it
#
# @param y the input GRanges
y_nclist_order <- function(y) {
  code <- as.integer(seqnames(y))
  y_start <- start(y)
  y_end <- end(y)
  if (anyNA(code) || any(y_end < y_start)) {
    return(FALSE)
  }
  if (is.unsorted(code)) {
    o <- order(code, method = "radix")
    code <- code[o]
    y_start <- y_start[o]
    y_end <- y_end[o]
  }
  n <- length(code)
  if (n < 2L) {
    return(TRUE)
  }
  k <- seq.int(2L, n)
  same <- code[k] == code[k - 1L]
  !any(same & (y_start[k] < y_start[k - 1L] |
    (y_start[k] == y_start[k - 1L] & y_end[k] > y_end[k - 1L])))
}

# The most bait blocks bootRanges() draws and moves at once, unless one
# replicate alone has more
boot_block_budget <- function() {
  BOOT_BLOCK_BUDGET
}

BOOT_BLOCK_BUDGET <- 1e6

# The number of bait blocks each replicate draws on the bait routes, as
# seg_bootstrap_prop(), seg_bootstrap_noprop() (before its blocks are cut
# to the tiling's) and unseg_across_draws() count them; NA where it cannot
# be computed. No draws.
#
# @param seg,proportionLength as in bootRanges()
# @param seg_mcols mcols(seg)
# @param chr_lens seqlengths(y)
# @param L_b the block length
# @param unseg_setup unseg_bootstrap_setup()'s value, without 'seg'
bait_blocks_per_rep <- function(seg, seg_mcols, chr_lens, L_b,
                                proportionLength, unseg_setup) {
  n <- tryCatch(suppressWarnings({
    if (is.null(seg)) {
      sum(ceiling(as.numeric(unseg_setup$chrom_lens) / unseg_setup$L_b))
    } else if (identical(as.vector(class(seg_mcols)), "DFrame")) {
      w <- as.numeric(width(seg))
      if (isTRUE(proportionLength)) {
        state <- seg_mcols@listData[["state"]]
        L_g <- sum(as.numeric(chr_lens))
        sum(vapply(seq_len(max(state)), function(s) {
          w_s <- w[state == s]
          sum(ceiling(w_s / round(L_b * sum(w_s) / L_g)))
        }, 0))
      } else {
        # the blocks kept: the tiling's of the segments in a state
        state <- seg_mcols@listData[["state"]]
        sum(ceiling(w[state %in% seq_len(max(state))] / L_b))
      }
    }
  }), error = function(e) NULL)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n)) {
    return(NA_real_)
  }
  n
}

# The replicates of each chunk: as many as boot_block_budget() holds, at
# least one, in order; all in one chunk where the count per replicate is NA
#
# @param n_rep the number of replicates, an integer
# @param per_rep the bait blocks per replicate, or NA
boot_chunks <- function(n_rep, per_rep) {
  k <- if (is.na(per_rep) || per_rep <= 0) n_rep else
    max(1, floor(boot_block_budget() / per_rep))
  if (k >= n_rep) {
    return(list(n_rep))
  }
  k <- as.integer(k)
  sizes <- rep.int(k, n_rep %/% k)
  if (n_rep %% k != 0L) {
    sizes <- c(sizes, n_rep %% k)
  }
  as.list(sizes)
}

# The moved ranges and counts per replicate of the chunks, in order, as one
# draw_move() of all replicates gives them; the 'block' column, which
# bootRanges() drops, is dropped first, as the chunks may differ in it
#
# @param parts per chunk, list(ranges, lens)
bind_chunks <- function(parts) {
  lens <- unlist(lapply(parts, `[[`, "lens"), use.names = FALSE)
  ranges <- lapply(parts, function(p) {
    x <- p$ranges
    if ("block" %in% names(mcols(x))) {
      mcols(x)$block <- NULL
    }
    x
  })
  rm(parts)
  list(ranges = do.call(c, ranges), lens = lens)
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
#
# With 'fuse', bulk_setup()'s value, and blocks$moved, the rows
# moved_positions() gives for these blocks (C_seg_draws_moved made them
# with the draws), the steps after moved_positions() run on them.
move_blocks_bulk <- function(y, blocks, exclude, excludeOption,
                             fuse = NULL) {
  if (is.list(blocks$moved) && !is.null(fuse)) {
    return(bulk_rows(fuse$y, blocks$moved, fuse$ex, exclude, excludeOption,
                     length(blocks$bait_len)))
  }
  n_blocks <- length(blocks$target_start)
  seg_code <- NULL
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
    # each segment's code; each block's, seg_code[blocks$bait_seg], is
    # built below only where moved_positions() does not derive it
    seg_code <- match(blocks$chr_names, seqlevels(y))
    bait_code <- NULL
  }
  circular <- unname(Seqinfo::isCircular(y)) %in% TRUE
  # moved_positions() takes bait_blocks()' inputs and the blocks' segments
  # where the check below could return NULL only for an NA code, on which
  # it declines: plain numbers, an integer segment per block and no
  # circular sequence. Where it gives no list, the steps below run as for
  # the other inputs.
  raw <- !is.null(seg_code) && !any(circular) &&
    plain_number(bait_start) && plain_number(blocks$width) &&
    plain_number(blocks$target_start) && is.integer(blocks$bait_seg) &&
    !is.object(blocks$bait_seg) &&
    length(blocks$bait_seg) == length(bait_start)
  if (!raw && !is.null(seg_code)) {
    bait_code <- seg_code[blocks$bait_seg]
  }
  if (!raw && bait_declines(y, bait_start, bait_code, blocks, circular)) {
    return(NULL)
  }
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
  # moved_positions() and moved_ranges() take y only as a GRanges, and
  # read its slots where the accessors would run updateObject() on them
  if (identical(as.vector(class(y)), "GRanges")) {
    y <- updateObject(y, check = FALSE)
  }
  moved <- NULL
  if (raw) {
    moved <- moved_positions(
      y, list(start = bait_start, width = blocks$width,
              target = blocks$target_start, seg = blocks$bait_seg,
              seg_code = seg_code),
      NULL, n_rep, NULL, target_chr, ex
    )
    if (!is.list(moved)) {
      bait_code <- seg_code[blocks$bait_seg]
      if (bait_declines(y, bait_start, bait_code, blocks, circular)) {
        return(NULL)
      }
    }
  }
  if (!is.list(moved)) {
    bait <- if (bait_plain(bait_start, bait_code, blocks)) {
      .Call(C_bait_blocks, bait_start, blocks$width, blocks$target_start,
            as.integer(n_rep))
    }
    if (is.null(bait)) {
      bait_ranges <- IRanges(
        start = bait_start,
        width = rep_len(blocks$width, length(bait_start))
      )
      block_shift <- rep(blocks$target_start, n_rep) - start(bait_ranges)
      bait <- bait_ranges
    } else {
      bait_ranges <- NULL
      block_shift <- bait$shift
    }
    moved <- moved_positions(
      y, bait, bait_code, n_rep, block_shift, target_chr, ex
    )
    if (isFALSE(moved)) {
      ex <- NULL
      moved <- moved_positions(
        y, bait, bait_code, n_rep, block_shift, target_chr
      )
    }
  }
  if (is.null(moved)) {
    ex <- NULL
    if (is.null(bait_ranges)) {
      # as built above from bait_start and the widths
      bait_ranges <- IRanges(start = bait$start, width = bait$width)
      block_shift <- rep(blocks$target_start, n_rep) - bait$start
    }
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
  bulk_rows(y, moved, ex, exclude, excludeOption, n_rep)
}

# move_blocks_bulk() where a replicate has at least as many bait blocks as
# 'y' has ranges, so that some sequence holds at least as many of its blocks
# as ranges of 'y' and moved_positions() declines: its route by
# bait_overlaps() and shift_swap_trim_positions(), run on pieces of at most
# 'piece' consecutive blocks, whose rows, in block order, are those of one
# pass. Each block's hits are found as bait_overlaps() finds them for all
# blocks at once: on the side its replicate's count on its sequence picks
# (the replicate's, counted over all pieces), and on y's NCList, prebuilt,
# wherever IRanges would build it on 'y'. Returns FALSE where the inputs are
# not of that kind (move_blocks_bulk() then runs), NULL where
# move_blocks_bulk() returns NULL, else its value.
#
# @param y the input GRanges
# @param blocks the draws, laid out by seg_draws_compiled(), without 'moved'
# @param exclude,excludeOption as in bootRanges()
# @param piece the most blocks moved at once
move_blocks_pieces <- function(y, blocks, exclude, excludeOption, piece) {
  len <- blocks$bait_len
  n_blocks <- length(blocks$target_start)
  start_all <- blocks$bait_start
  seg_all <- blocks$bait_seg
  width <- blocks$width
  target_chr <- blocks$target_chr
  # a factor whose subsets keep its attributes, as shift_swap_trim_positions()
  # reads them
  chr_fac <- if (is(target_chr, "Rle")) runValue(target_chr) else target_chr
  if (!identical(as.vector(class(y)), "GRanges") || !is.null(blocks$moved) ||
      is.null(start_all) || !is.integer(seg_all) || is.object(seg_all) ||
      length(seg_all) != length(start_all) || length(len) < 1L ||
      !plain_number(start_all) || !plain_number(width) ||
      !(length(width) %in% c(1L, n_blocks)) ||
      !plain_number(blocks$target_start) || !is.factor(chr_fac) ||
      !setequal(names(attributes(chr_fac)), c("levels", "class")) ||
      any(unname(Seqinfo::isCircular(y)) %in% TRUE) ||
      !any(len >= length(y))) {
    return(FALSE)
  }
  n_rep <- length(len)
  if (any(len != n_blocks)) {
    return(NULL)
  }
  lev <- seqlevels(y)
  n_lev <- length(lev)
  seg_code <- match(blocks$chr_names, lev)
  n <- length(start_all)
  size <- as.integer(max(1, min(floor(piece), n)))
  bounds <- unique(c(seq.int(0L, n, by = size), n))
  piece_at <- function(k) seq.int(bounds[k] + 1L, bounds[k + 1L])
  rep_of <- function(i) (i - 1L) %/% n_blocks + 1L
  # each replicate's blocks per sequence, and bait_declines() over all
  n_bait <- numeric(n_rep * n_lev)
  for (k in seq_len(length(bounds) - 1L)) {
    i <- piece_at(k)
    code <- seg_code[seg_all[i]]
    if (anyNA(code)) {
      return(NULL)
    }
    n_bait <- n_bait + tabulate((rep_of(i) - 1L) * n_lev + code, n_rep * n_lev)
  }
  y <- updateObject(y, check = FALSE)
  n_y <- tabulate(as.integer(seqnames(y)), n_lev)
  y_nclist <- NULL
  parts <- vector("list", length(bounds) - 1L)
  n_hit <- 0
  for (k in seq_along(parts)) {
    i <- piece_at(k)
    t <- (i - 1L) %% n_blocks + 1L
    code <- seg_code[seg_all[i]]
    bait_ranges <- IRanges(
      start = start_all[i],
      width = if (length(width) == 1L) rep_len(width, length(i)) else width[t]
    )
    block_shift <- blocks$target_start[t] - start(bait_ranges)
    random_blocks <- GRanges(seqnames = blocks$chr_names[seg_all[i]],
                             ranges = bait_ranges)
    rep_chr <- (rep_of(i) - 1L) * n_lev + code
    on_y <- n_bait[rep_chr] >= n_y[code]
    if (any(on_y) && is.null(y_nclist)) {
      y_nclist <- GNCList(y)
    }
    fo <- piece_overlaps(random_blocks, code, which(!on_y), which(on_y), y,
                         y_nclist)
    hit_block <- fo$query
    tc <- target_chr[t]
    if (is(tc, "Rle")) {
      tc <- decode(tc)
    }
    n_hit <- n_hit + length(hit_block)
    if (length(hit_block) == 0L) {
      next
    }
    moved <- shift_swap_trim_positions(
      y, fo$subject, block_shift[hit_block], tc, 1L, hit_block
    )
    if (is.null(moved)) {
      return(NULL)
    }
    w <- which(moved$width > 0L)
    parts[[k]] <- list(
      s = moved$s[w], code = moved$code[w], start = moved$start[w],
      width = moved$width[w], block = rep_of(i[hit_block[w]]),
      chr_attr = moved$chr_attr
    )
  }
  if (n_hit == 0) {
    return(NULL)
  }
  parts <- parts[!vapply(parts, is.null, NA)]
  moved <- list(
    s = unlist(lapply(parts, `[[`, "s"), use.names = FALSE),
    code = unlist(lapply(parts, `[[`, "code"), use.names = FALSE),
    start = unlist(lapply(parts, `[[`, "start"), use.names = FALSE),
    width = unlist(lapply(parts, `[[`, "width"), use.names = FALSE),
    block = unlist(lapply(parts, `[[`, "block"), use.names = FALSE),
    chr_attr = parts[[1L]]$chr_attr
  )
  rm(parts)
  bulk_rows(y, moved, NULL, exclude, excludeOption, n_rep)
}

# bait_overlaps() for one piece of the blocks, given which of them find
# their hits on the bait side ('on_bait') and which on y's NCList ('on_y',
# prebuilt as 'y_nclist'), by block and within a block as bait_overlaps()
# orders them
#
# @param bait the piece's bait blocks
# @param bait_code the index of each block's sequence in seqlevels(y)
# @param on_bait,on_y the blocks of each kind, by index in 'bait'
# @param y the input GRanges
# @param y_nclist GNCList(y), or NULL when 'on_y' is empty
piece_overlaps <- function(bait, bait_code, on_bait, on_y, y, y_nclist) {
  query <- integer()
  subject <- integer()
  if (length(on_bait) != 0L) {
    fo <- sorted_overlaps(bait, bait_code, on_bait, y)
    if (is.null(fo)) {
      fo <- findOverlaps(GNCList(bait[on_bait]), y)
      fo <- list(query = on_bait[queryHits(fo)], subject = subjectHits(fo))
    }
    query <- fo$query
    subject <- fo$subject
  }
  if (length(on_y) != 0L) {
    fo_y <- findOverlaps(bait[on_y], y_nclist)
    query <- c(query, on_y[queryHits(fo_y)])
    subject <- c(subject, subjectHits(fo_y))
    if (length(on_bait) != 0L) {
      o <- order(query, method = "radix")
      query <- query[o]
      subject <- subject[o]
    }
  }
  list(query = query, subject = subject)
}

# move_blocks_bulk()'s steps after the moved positions: the exclude step,
# if moved_positions() did not apply it, and the ranges with the number
# per replicate
#
# @param y the input GRanges, updated as move_blocks_bulk() updates it
# @param moved moved_positions()' list, or the list built from the hits
# @param ex the exclude tables moved_positions() was given, or NULL
# @param exclude,excludeOption as in bootRanges()
# @param n_rep the number of replicates
bulk_rows <- function(y, moved, ex, exclude, excludeOption, n_rep) {
  # only rows of positive width from here: overlapsAny() is per range
  s <- moved$s
  # NULL when moved_positions() got 'ex': moved_ranges() reads its runs
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

# What move_blocks_bulk() computes from 'y' and 'exclude' alone on its way
# to moved_positions() with bait_blocks()' inputs (its 'raw' route): the
# updated 'y', the exclude tables 'ex' (NULL where exclude_tables() fails
# or warns), the seqlevels of 'y' before and after the update ('lev0',
# 'lev') and moved_positions()' arguments to C_moved_rows_bait_par that do
# not depend on the blocks ('args', and 'threads', omp_threads() read once
# for every chunk's compiled draws). NULL where that route is not taken
# whatever the blocks: 'y' not a GRanges, a circular sequence, too many
# seqlevels, or seqlevels the update changes.
#
# @param y the input GRanges
# @param exclude,excludeOption as in bootRanges()
bulk_setup <- function(y, exclude, excludeOption) {
  if (!identical(as.vector(class(y)), "GRanges") ||
      any(unname(Seqinfo::isCircular(y)) %in% TRUE)) {
    return(NULL)
  }
  if (!is.null(exclude)) {
    ex <- tryCatch(
      exclude_tables(y, exclude, excludeOption),
      warning = function(w) NULL, error = function(e) NULL
    )
  } else {
    ex <- list(mode = 0L, seg = integer(), start = integer(),
               end = integer(), seqinfo = NULL)
  }
  lev0 <- seqlevels(y)
  y <- updateObject(y, check = FALSE)
  lev <- seqlevels(y)
  if (length(lev) >= 2097152L || !identical(lev, lev0)) {
    return(NULL)
  }
  circular <- unname(Seqinfo::isCircular(y)) %in% TRUE
  seqlen <- unname(seqlengths(y))
  cols <- NULL
  ex_c <- ex
  with_code <- is.null(ex)
  if (with_code) {
    ex_c <- list(mode = 0L, seg = integer(), start = integer(),
                 end = integer())
  } else {
    cols <- gather_columns(y)
  }
  rng <- y@ranges
  y_start <- start(rng)
  y_width <- width(rng)
  list(
    y = y, ex = ex, lev0 = lev0, lev = lev, threads = omp_threads(),
    args = list(
      as.integer(seqnames(y)), y_start, y_width - 1L + y_start, y_width,
      circular, as.integer(seqlen), !(circular | is.na(seqlen)),
      ex_c$mode, ex_c$seg, ex_c$start, ex_c$end, cols, with_code
    )
  )
}

# The arguments C_seg_draws_moved_par takes besides the draws' for blocks of
# these widths, targets and segment sequences, as move_blocks_bulk() and
# moved_positions() would build them; NULL where move_blocks_bulk() would
# not take its 'raw' route to moved_positions() for any draws.
#
# @param fuse bulk_setup()'s value
# @param chr_names the sequence of each segment the draws index
# @param width the block widths
# @param target_start,target_chr the blocks' targets, as the draws' value
# holds them (target_chr a factor or an Rle of one)
fuse_args <- function(fuse, chr_names, width, target_start, target_chr) {
  if (is(target_chr, "Rle")) {
    target_chr <- decode(target_chr)
  }
  chr_attr <- attributes(target_chr)
  if (!plain_number(width) || !plain_number(target_start) ||
      !is.factor(target_chr) ||
      !setequal(names(chr_attr), c("levels", "class")) ||
      !identical(levels(target_chr), fuse$lev)) {
    return(NULL)
  }
  list(
    seg_code = match(chr_names, fuse$lev0), width = width,
    target_start = target_start, t_code = as.integer(target_chr),
    chr_attr = chr_attr, args = fuse$args, threads = fuse$threads
  )
}

# Plain numbers: then bait_declines() is FALSE without a circular
# sequence, and bait_blocks() gives what moved_positions() reads of the
# IRanges and the shifts, wherever building them could not fail or warn
#
# @param bait_start,bait_code,blocks as in move_blocks_bulk()
bait_plain <- function(bait_start, bait_code, blocks) {
  plain_number(bait_start) && plain_number(blocks$width) &&
    plain_number(blocks$target_start) &&
    length(bait_code) == length(bait_start)
}

# TRUE where move_blocks_bulk() returns NULL before moving any block: a
# block on a sequence 'y' lacks, or starting outside a circular one
#
# @param y the input GRanges
# @param bait_start,bait_code,blocks as in move_blocks_bulk()
# @param circular isCircular(y), FALSE for NA
bait_declines <- function(y, bait_start, bait_code, blocks, circular) {
  anyNA(bait_code) ||
    ((!bait_plain(bait_start, bait_code, blocks) || any(circular)) &&
      any(circular[bait_code] & (bait_start < 1 |
      bait_start > unname(seqlengths(y))[bait_code])))
}

# An integer or double vector without a class
#
# @param x any object
plain_number <- function(x) {
  (is.integer(x) || is.double(x)) && !is.object(x)
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
# of 's', when gather_columns() gives them, and 'code' is NULL: only its
# runs are returned, which moved_ranges() reads in its place.
#
# @param y the input GRanges, updated (updateObject()), whose slots are read
# @param bait the ranges (IRanges) of the bait blocks of all replicates,
# in replicate order, or bait_blocks()' list of their integer starts and
# widths, or bait_blocks()' inputs ('start', 'width', 'target') with each
# block's segment ('seg') and each segment's index in seqlevels(y)
# ('seg_code'), from which the C routine derives the blocks' starts,
# widths, shifts and codes, declining wherever bait_blocks() would or a
# code is NA
# @param bait_code the index of each block's sequence in seqlevels(y);
# unread with 'seg_code'
# @param n_rep the number of replicates
# @param block_shift the shift of each block; unread with 'seg_code'
# @param target_chr the sequence each block moves to, a factor
# @param ex the output of exclude_tables(), the tables of no exclude step
# (mode 0), or NULL
moved_positions <- function(y, bait, bait_code, n_rep, block_shift,
                            target_chr, ex = NULL) {
  chr_attr <- attributes(target_chr)
  lev <- seqlevels(y)
  raw <- is.list(bait) && !is.null(bait$seg_code)
  if (!identical(as.vector(class(y)), "GRanges") ||
      (!raw && !is.numeric(block_shift)) || length(lev) >= 2097152L ||
      !is.factor(target_chr) ||
      !setequal(names(chr_attr), c("levels", "class")) ||
      !identical(levels(target_chr), lev)) {
    return(NULL)
  }
  circular <- unname(Seqinfo::isCircular(y)) %in% TRUE
  seqlen <- unname(seqlengths(y))
  cols <- NULL
  with_code <- is.null(ex)
  if (with_code) {
    ex <- list(mode = 0L, seg = integer(), start = integer(), end = integer())
  } else {
    cols <- gather_columns(y)
  }
  # move_blocks_bulk() updated y, so its ranges slot is ranges(y); end()
  # computed as the Ranges method does
  rng <- y@ranges
  y_start <- start(rng)
  y_width <- width(rng)
  moved <- if (raw) {
    .Call(
      C_moved_rows_bait_par, bait$seg, bait$seg_code, bait$start,
      bait$width, as.integer(n_rep), as.integer(seqnames(y)), y_start,
      y_width - 1L + y_start, y_width,
      circular, as.integer(seqlen), !(circular | is.na(seqlen)),
      bait$target, as.integer(target_chr),
      ex$mode, ex$seg, ex$start, ex$end, cols, with_code, omp_threads()
    )
  } else {
    if (is.list(bait)) {
      b_start <- bait$start
      b_width <- bait$width
    } else {
      b_start <- start(bait)
      b_width <- width(bait)
    }
    .Call(
      C_moved_rows_opt_par, as.integer(bait_code), b_start, b_width,
      as.integer(n_rep), as.integer(seqnames(y)), y_start,
      y_width - 1L + y_start, y_width,
      circular, as.integer(seqlen), !(circular | is.na(seqlen)),
      suppressWarnings(as.integer(block_shift)), as.integer(target_chr),
      ex$mode, ex$seg, ex$start, ex$end, cols, with_code, omp_threads()
    )
  }
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
# @param y the input GRanges, updated (updateObject()), whose slots are read
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
      y@ranges, start = start, width = width, check = FALSE
    )
    ans_mcols <- y@elementMetadata
    if (is.null(cols)) {
      cols <- lapply(as.list(ans_mcols), `[`, s)
    }
    ans_mcols <- BiocGenerics:::replaceSlots(
      ans_mcols, listData = cols, nrows = n, check = FALSE
    )
  } else {
    ans_ranges <- update_ranges(
      extractROWS(y@ranges, s), start = start, width = width
    )
    ans_mcols <- extractROWS(y@elementMetadata, s)
  }
  if (!is.null(block)) {
    ans_mcols$block <- block
  }
  # update() without its updateObject(y), and with the seqinfo in the
  # same call
  BiocGenerics:::replaceSlots(
    y, seqnames = ans_seqnames, ranges = ans_ranges, strand = ans_strand,
    elementMetadata = ans_mcols,
    .slotList = if (!is.null(seqinfo)) list(seqinfo = seqinfo),
    check = FALSE
  )
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
  si_y <- seqinfo(y)
  si_ex <- seqinfo(exclude)
  # merge() of two identical Seqinfo objects warns of nothing and keeps
  # their circularity
  si <- if (identical(si_y, si_ex)) si_y else merge(si_y, si_ex)
  lev <- seqlevels(y)
  # one ranges() call for start(), width() and end() (the Ranges method)
  rng <- ranges(exclude)
  rng_start <- start(rng)
  rng_width <- width(rng)
  if (length(lev) >= 2097152L || any(Seqinfo::isCircular(si) %in% TRUE) ||
      !all(runValue(strand(exclude)) == "*") || any(rng_width == 0L)) {
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
    ex_start <- key + rng_start[on_y]
    o <- order(ex_start, method = "radix")
    ex_start <- ex_start[o]
    ex_end <- cummax((key[o] + (rng_width - 1L + rng_start)[on_y][o]))
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
    tb <- interval_gaps(y, exclude)
    if (is.null(tb)) {
      gap <- gaps(exclude, end = seqlengths(y))
      gap <- gap[strand(gap) == "*"]
      tb <- gap_tables(y, gap)
      if (is.null(tb) || !tb$ordered) {
        return(NULL)
      }
    }
    return(list(mode = 2L, seg = tb$seg, start = tb$start, end = tb$end,
                seqinfo = tb$seqinfo))
  }
  NULL
}

# gap_tables(y, gap) for gap <- gaps(exclude, end = seqlengths(y)) on
# strand "*", built from exclude_intervals(): per sequence of 'y', the parts
# of [1, seqlength] between its merged intervals (before the first, between
# two, after the last; all of it when it has none), those of positive
# width, in order. NULL unless 'y' and 'exclude' are plain GRanges with
# identical seqinfo, every seqlength is known and below the integer
# maximum, exclude_intervals() takes 'exclude' (all on "*", no empty range,
# no circular sequence) and every range of 'exclude' lies within [1, its
# seqlength]; gaps() then makes the same gaps, and gap_tables() keeps
# seqinfo(y).
#
# @param y the input GRanges
# @param exclude the excluded regions
interval_gaps <- function(y, exclude) {
  if (!identical(as.vector(class(y)), "GRanges") ||
      !identical(as.vector(class(exclude)), "GRanges") ||
      !identical(seqinfo(exclude), seqinfo(y))) {
    return(NULL)
  }
  seqlen <- unname(seqlengths(y))
  # IRanges' gaps() with end = .Machine$integer.max ignores the ranges
  # (IRanges 2.47.5: [1, end] plus a gap after each range)
  if (anyNA(seqlen) || any(seqlen == .Machine$integer.max)) {
    return(NULL)
  }
  iv <- exclude_intervals(y, exclude)
  if (is.null(iv)) {
    return(NULL)
  }
  n_lev <- length(seqlen)
  n_int <- diff(iv$seg)
  int_code <- rep.int(seq_len(n_lev), n_int)
  if (any(iv$start < 1L) || any(iv$end > seqlen[int_code])) {
    return(NULL)
  }
  # sequence c holds candidates seg[c] + c to seg[c + 1] + c: the one
  # before interval i is i + c - 1, the one after it i + c (doubles, so
  # that an end at the integer maximum does not overflow)
  n_cand <- length(int_code) + n_lev
  cand_code <- rep.int(seq_len(n_lev), n_int + 1L)
  cand_start <- rep.int(1, n_cand)
  cand_end <- as.double(seqlen)[cand_code]
  at <- seq_along(int_code) + int_code
  cand_start[at] <- as.double(iv$end) + 1
  cand_end[at - 1L] <- as.double(iv$start) - 1
  keep <- which(cand_start <= cand_end)
  gap_code <- cand_code[keep]
  list(
    seg = c(0L, cumsum(tabulate(gap_code, n_lev))),
    start = as.integer(cand_start[keep]), end = as.integer(cand_end[keep]),
    seqinfo = seqinfo(y)
  )
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
  rng <- ranges(gap)
  gap_start <- start(rng)[on_x][o]
  gap_end <- end(rng)[on_x][o]
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
                                proportion_length = TRUE, fuse = NULL) {

  # sequence along the states: 1,2,..,S
  esses <- seq_len(max(state))
  chr_names <- as.character(chr_names)

  # break apart this into two sub-functions, defined below
  if (proportion_length) {
    out <- seg_bootstrap_prop(seg, state, L_b, chr_names, chr_lens,
                              esses, n_rep, fuse)
  } else {
    out <- seg_bootstrap_noprop(seg, state, L_b, chr_names, chr_lens,
                                esses, n_rep, fuse)
  }
  # with 'fuse' the compiled branches made it before the draws
  if (!isTRUE(out$chr_factor)) {
    out$target_chr <- factor(out$target_chr, levels = seqlevels(y))
  }
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
# With 'fz', fuse_args()' value, C_seg_draws_moved makes the same draws and
# the rows moved_positions() gives for them, as a fourth element.
#
# @param tabs per table, list(index, width, start, nblocks, L, n)
# @param n_rep the number of replicates
# @param st the state of each segment (0 for none), noprop only
# @param rearr the number of tiling blocks per state, noprop only
# @param fz fuse_args()' value, or NULL
seg_draws_compiled <- function(tabs, n_rep, st = NULL, rearr = NULL,
                               fz = NULL) {
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
  prob <- as.double(unlist(lapply(tabs, `[[`, "width"), use.names = FALSE))
  rmin <- as.double(unlist(lapply(tabs, `[[`, "start"), use.names = FALSE))
  rmax <- unlist(rmax, use.names = FALSE)
  gidx <- as.integer(unlist(lapply(tabs, `[[`, "index"), use.names = FALSE))
  off <- c(0L, cumsum(m))
  if (is.null(fz)) {
    return(.Call(C_seg_draws, n_rep, prob, rmin, rmax, gidx, off, n_draw,
                 st, rearr))
  }
  a <- fz$args
  .Call(
    C_seg_draws_moved_par, n_rep, prob, rmin, rmax, gidx, off, n_draw, st,
    rearr, fz$seg_code, fz$width, a[[1L]], a[[2L]], a[[3L]], a[[4L]],
    a[[5L]], a[[6L]], a[[7L]], fz$target_start, fz$t_code, a[[8L]],
    a[[9L]], a[[10L]], a[[11L]], a[[12L]], a[[13L]], fz$threads
  )
}

# The draws' value with the moved rows of C_seg_draws_moved, if any, as
# moved_positions() returns them
#
# @param d seg_draws_compiled()'s value
# @param fz the fuse_args() it was given, or NULL
draws_moved <- function(d, fz) {
  moved <- if (length(d) > 3L) d[[4L]]
  if (!is.list(moved)) {
    return(NULL)
  }
  moved$chr_attr <- fz$chr_attr
  moved
}

# breaking up the seg_bootstrap_draws() code... proportional length case
seg_bootstrap_prop <- function(seg, state, L_b, chr_names, chr_lens,
                               esses, n_rep, fuse = NULL) {
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

  out <- seg_prop_compiled(per_state, L_b, L_b_per_state, n_rep, fuse,
                           chr_names)
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
# computed first; NULL where the R loop must run. Where fuse_args() takes
# 'fuse' (with chr_names, the segments' sequences), target_chr is made as
# seg_bootstrap_draws() makes it, before the draws, and the draws move the
# blocks (fused_parts()).
seg_prop_compiled <- function(per_state, L_b, L_b_per_state, n_rep,
                              fuse = NULL, chr_names = NULL) {
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
  target_chr <- NULL
  fz <- NULL
  if (!is.null(fuse)) {
    fz <- tryCatch({
      target_chr <- factor(do.call(c, lapply(tiles, `[[`, 2)),
                           levels = fuse$lev0)
      nb <- vapply(tabs, function(tb) as.integer(tb$n), 1L)
      fuse_args(fuse, chr_names, rep(L_b_per_state, nb),
                unlist(lapply(tiles, `[[`, 1)), target_chr)
    }, warning = function(w) NULL, error = function(e) NULL)
  }
  d <- seg_draws_compiled(tabs, n_rep, fz = fz)
  if (is.null(d)) {
    return(NULL)
  }
  nblocks_per_state <- vapply(tabs, function(tb) as.integer(tb$n), 1L)
  out <- list(
    draws = NULL,
    width = rep(L_b_per_state, nblocks_per_state),
    target_start = unlist(lapply(tiles, `[[`, 1)),
    target_chr = if (is.null(fz)) do.call(c, lapply(tiles, `[[`, 2)) else
      target_chr,
    bait_start = d[[1L]], bait_seg = d[[2L]], bait_len = d[[3L]]
  )
  fused_parts(out, d, fz)
}

# 'out' with, after C_seg_draws_moved, 'chr_factor' (TRUE: target_chr is
# already the factor seg_bootstrap_draws() makes) and 'moved', the rows of
# draws_moved() where it gives any; 'out' itself without 'fz'
#
# @param out the draws' value
# @param d,fz as in draws_moved()
fused_parts <- function(out, d, fz) {
  if (!is.null(fz)) {
    out$chr_factor <- TRUE
    out$moved <- draws_moved(d, fz)
  }
  out
}

# breaking up the seg_bootstrap_draws() code... non-proportional length case
seg_bootstrap_noprop <- function(seg, state, L_b, chr_names, chr_lens,
                                 esses, n_rep, fuse = NULL) {

  seg_width <- width(seg)
  seg_start <- start(seg)
  # number of rearranged blocks within each range of the segmentation
  nblocks_per_seg <- ceiling(seg_width / L_b)
  # total number of tiling blocks
  n <- sum(nblocks_per_seg)

  out <- seg_noprop_compiled(seg_width, seg_start, nblocks_per_seg, n,
                             L_b, state, chr_names, esses, n_rep, fuse)
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
# computed first; NULL where the R loop must run. 'fuse' as in
# seg_prop_compiled().
seg_noprop_compiled <- function(seg_width, seg_start, nblocks_per_seg, n,
                                L_b, state, chr_names, esses, n_rep,
                                fuse = NULL) {
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
  target_chr <- NULL
  fz <- NULL
  if (!is.null(fuse)) {
    fz <- tryCatch({
      target_chr <- factor(do.call(c, tiles[[2L]]), levels = fuse$lev0)
      fuse_args(fuse, chr_names, L_b, unlist(tiles[[1L]]), target_chr)
    }, warning = function(w) NULL, error = function(e) NULL)
  }
  d <- seg_draws_compiled(tabs, n_rep,
                          st = match(state, esses, nomatch = 0L),
                          rearr = lengths(tiles[[1L]]), fz = fz)
  if (is.null(d)) {
    return(NULL)
  }
  out <- list(
    draws = NULL,
    width = L_b,
    target_start = unlist(tiles[[1L]]),
    target_chr = if (is.null(fz)) do.call(c, tiles[[2L]]) else target_chr,
    bait_start = d[[1L]], bait_seg = d[[2L]], bait_len = d[[3L]]
  )
  fused_parts(out, d, fz)
}
