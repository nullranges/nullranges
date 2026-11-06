library(nullranges)
library(GenomicRanges)

## The compiled moves give each bait block's hits by increasing index in
## 'y'. findOverlaps() gives them in that order where IRanges builds its
## NCList on the blocks, and in the NCList's order (start ascending, end
## descending) where a replicate has at least as many blocks on a sequence
## as 'y' has ranges there. So the compiled moves take a sequence with no
## ranges, or whose ranges are already in NCList order, and leave to the R
## route only a sequence holding an equal-start pair whose later range ends
## later. Every call is compared with base's code under the same seed.

## the value (or error message), warnings and .Random.seed after f() from
## set.seed(seed)
ncOutcome <- function(f, seed = 1L) {
  set.seed(seed)
  warns <- character()
  value <- tryCatch(withCallingHandlers(f(), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) conditionMessage(e))
  list(value = value, warnings = warns,
       seed = sprintf("%d", get(".Random.seed", envir = globalenv())))
}

## the outcome of bootRanges() and the routes it took: "compiled" where the
## compiled moves gave the rows (with the draws or after them), "R" where
## the blocks were moved by findOverlaps() (bait_overlaps() or
## move_blocks_pieces())
ncRouted <- function(args) {
  ns <- asNamespace("nullranges")
  dm <- get("draws_moved", ns)
  mp <- get("moved_positions", ns)
  bo <- get("bait_overlaps", ns)
  mbp <- get("move_blocks_pieces", ns)
  route <- character()
  local_mocked_bindings(
    draws_moved = function(...) {
      r <- dm(...)
      if (!is.null(r)) route <<- c(route, "compiled")
      r
    },
    moved_positions = function(...) {
      r <- mp(...)
      if (is.list(r)) route <<- c(route, "compiled")
      r
    },
    bait_overlaps = function(...) {
      route <<- c(route, "R")
      bo(...)
    },
    move_blocks_pieces = function(...) {
      r <- mbp(...)
      if (is.list(r)) route <<- c(route, "R")
      r
    },
    .package = "nullranges"
  )
  out <- ncOutcome(function() do.call(bootRanges, args))
  out$route <- unique(route)
  out
}

## 'n1' sorted ranges on chr1 (1e7), 'y2' on chr2 (1e6) and a chrY of 1e6
## with none; a two-state segmentation of chr1 and chr2 (none of chrY) and
## an exclusion list with intervals on chr1, on the empty chrY and on chrZ,
## a sequence 'y' lacks
ncInputs <- function(y2, n1 = 2000L) {
  set.seed(11)
  len <- c(chr1 = 1e7, chr2 = 1e6, chrY = 1e6)
  y1 <- IRanges(sort(sample(1e7 - 5000, n1)),
                width = sample(50:2000, n1, replace = TRUE))
  y <- GRanges(c(rep("chr1", n1), rep("chr2", length(y2))), c(y1, y2),
               seqlengths = len)
  y$score <- seq_along(y)
  seg <- GRanges(c("chr1", "chr1", "chr2"),
                 IRanges(c(1, 5e6 + 1, 1), c(5e6, 1e7, 1e6)),
                 state = c(1L, 2L, 1L))
  seqlevels(seg) <- names(len)
  seqlengths(seg) <- len
  ex <- GRanges(c("chr1", "chrY", "chrZ"),
                IRanges(c(2e6, 1e5, 1e3), width = c(1e5, 5e4, 1e4)),
                seqlengths = c(len, chrZ = 1e5))
  list(y = y, seg = seg, ex = ex)
}

ncCalls <- function(inp) {
  list(
    across = list(),
    drop = list(exclude = inp$ex),
    seg_drop = list(seg = inp$seg, exclude = inp$ex),
    noprop_trim = list(seg = inp$seg, proportionLength = FALSE,
                       exclude = inp$ex, excludeOption = "trim")
  )
}

## each call against base; 'expect_route' the route expected where base
## returns a value, or NULL where only the outcome is checked
ncCheck <- function(inp, L, R, label, expect_route = NULL) {
  calls <- ncCalls(inp)
  for (nm in names(calls)) {
    args <- c(list(y = inp$y, blockLength = L, R = R), calls[[nm]])
    base <- ncOutcome(function() {
      do.call(nullranges:::bootRanges_base, args)
    })
    got <- ncRouted(args)
    info <- paste(label, nm)
    expect_identical(got$value, base$value, info = info)
    expect_identical(got$warnings, base$warnings, info = info)
    expect_identical(got$seed, base$seed, info = info)
    if (is.character(base$value)) {
      next
    }
    if (identical(expect_route, "compiled")) {
      expect_true("compiled" %in% got$route, info = info)
      expect_false("R" %in% got$route, info = info)
    } else if (identical(expect_route, "R")) {
      expect_true("R" %in% got$route, info = info)
    }
  }
}

test_that("an empty sequence keeps the compiled moves", {
  inp <- ncInputs(IRanges())
  ncCheck(inp, 1e5, 20, "empty chr2", "compiled")
})

test_that("a sequence with fewer ranges than blocks in NCList order", {
  ## distinct starts, and an equal-start pair whose later range is shorter
  ## (not sorted, so the calls without 'seg' give base's error)
  y2 <- IRanges(c(1000, 200000, 200000, 450000, 451000),
                c(80000, 260000, 210000, 470000, 452000))
  inp <- ncInputs(y2)
  ncCheck(inp, 1e5, 20, "ordered chr2", "compiled")
})

test_that("an equal-start pair whose later range ends later", {
  ## sorted GRanges put the shorter of an equal-start pair first; y's NCList
  ## visits the longer first, so these blocks keep the R route
  y2 <- IRanges(c(1000, 200000, 200000, 450000, 451000),
                c(80000, 210000, 260000, 470000, 452000))
  inp <- ncInputs(y2)
  expect_false(nullranges:::y_nclist_order(inp$y))
  ncCheck(inp, 1e5, 20, "equal start chr2", "R")
})

test_that("many short contigs with a few ranges each", {
  set.seed(21)
  n_ctg <- 300L
  len <- setNames(sample(5e4:2e5, n_ctg, replace = TRUE),
                  sprintf("ctg%03d", seq_len(n_ctg)))
  k <- sample(0:3, n_ctg, replace = TRUE)
  chr <- rep(names(len), k)
  st <- unlist(lapply(seq_len(n_ctg), function(i) {
    sort(sample(len[[i]] - 3000, k[i]))
  }))
  y <- GRanges(chr, IRanges(st, width = sample(100:2000, length(st), TRUE)),
               seqlengths = len)
  seqlevels(y) <- names(len)
  seqlengths(y) <- len
  y$score <- seq_along(y)
  half <- names(len)[seq_len(n_ctg / 2)]
  seg <- GRanges(names(len), IRanges(1, len),
                 state = rep(1:2, length.out = n_ctg), seqlengths = len)
  ex <- GRanges(half[1:3], IRanges(1000, width = 2000), seqlengths = len)
  inp <- list(y = y, seg = seg, ex = ex)
  expect_true(nullranges:::y_nclist_order(y))
  ncCheck(inp, 5e4, 10, "contigs", "compiled")
})

test_that("a small y with fewer ranges than blocks per replicate", {
  set.seed(31)
  len <- c(chr1 = 1e6, chr2 = 8e5, chrY = 5e5)
  y <- GRanges(rep(c("chr1", "chr2"), c(20, 10)),
               IRanges(c(sort(sample(9e5, 20)), sort(sample(7e5, 10))),
                       width = sample(50:5000, 30, replace = TRUE)),
               seqlengths = len)
  y$score <- seq_along(y)
  seg <- GRanges(c("chr1", "chr1", "chr2"),
                 IRanges(c(1, 500001, 1), c(5e5, 1e6, 8e5)),
                 state = c(1L, 2L, 1L))
  seqlevels(seg) <- names(len)
  seqlengths(seg) <- len
  ex <- GRanges(c("chr1", "chrY"), IRanges(c(1e5, 1e4), width = 5e4),
                seqlengths = len)
  inp <- list(y = y, seg = seg, ex = ex)
  ncCheck(inp, 1e4, 5, "small y", "compiled")

  ## the same with an equal-start pair whose later range ends later: the
  ## draws are made without the compiled moves, which would decline them
  long <- y[c(3, 25)]
  end(long) <- end(long) + 3000
  inp$y <- sort(c(y, long))
  expect_false(nullranges:::y_nclist_order(inp$y))
  ncCheck(inp, 1e4, 5, "small y, equal start", "R")
})

test_that("y_nclist_order reads each sequence in index order", {
  y <- GRanges(c("chr2", "chr1", "chr2", "chr1"),
               IRanges(c(5, 10, 5, 10), c(20, 30, 10, 40)))
  ## chr1: (10, 30) then (10, 40): the later one ends later
  expect_false(nullranges:::y_nclist_order(y))
  expect_true(nullranges:::y_nclist_order(y[c(1, 2, 3)]))
  expect_true(nullranges:::y_nclist_order(y[c(4, 2, 1, 3)]))
  ## starts out of order within a sequence
  z <- GRanges("chr1", IRanges(c(20, 10), width = 5))
  expect_false(nullranges:::y_nclist_order(z))
  expect_true(nullranges:::y_nclist_order(GRanges()))
})
