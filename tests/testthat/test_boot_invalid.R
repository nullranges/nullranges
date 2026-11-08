library(nullranges)
library(GenomicRanges)

## bootRanges() on input it rejects, or rejects in some replicate. Each
## replicate draws its bait blocks, then moves them, before the next draws,
## and an error stops the loop, so the error comes at a fixed point of the
## random stream: these tests require the error, its warnings and the
## .Random.seed left behind, derived from that order.

## the error (or NULL), the value, the warnings and .Random.seed after f()
## from set.seed(seed)
bootOutcome <- function(f, seed = 1L) {
  set.seed(seed)
  warns <- character()
  value <- NULL
  err <- tryCatch(withCallingHandlers({
    value <- f()
    NULL
  }, warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(error = err, value = value, warnings = warns,
       seed = get(".Random.seed", envir = globalenv()))
}

## two .Random.seed vectors identical (compared as text: testthat's numeric
## comparison overflows on their differences when reporting a mismatch)
expectSameSeed <- function(object, expected, info = NULL) {
  expect_identical(sprintf("%d", object), sprintf("%d", expected),
                   info = info)
}

## .Random.seed after set.seed(seed) and nothing else
seedAfterSet <- function(seed = 1L) {
  set.seed(seed)
  get(".Random.seed", envir = globalenv())
}

## 80 ranges on two sequences, a two-state segmentation and an excluded
## region
invalidInputs <- function() {
  set.seed(3)
  y <- GRanges(rep(c("chr1", "chr2"), each = 40),
               IRanges(c(sort(sample(1e6, 40)), sort(sample(8e5, 40))),
                       width = 200),
               seqlengths = c(chr1 = 1e6 + 500, chr2 = 8e5 + 500))
  y$score <- seq_along(y)
  seg <- GRanges(c("chr1", "chr1", "chr2"),
                 IRanges(c(1, 500001, 1), c(5e5, 1e6 + 500, 8e5 + 500)),
                 state = c(1L, 2L, 1L), seqlengths = seqlengths(y))
  ex <- GRanges("chr1", IRanges(1e5, 2e5), seqlengths = seqlengths(y))
  list(y = y, seg = seg, ex = ex)
}

## every way bootRanges() moves blocks, by its arguments
bootRoutes <- function(seg) {
  list(across = list(),
       permute = list(type = "permute"),
       within = list(withinChrom = TRUE),
       permuteWithin = list(type = "permute", withinChrom = TRUE),
       segmented = list(seg = seg),
       segmentedFixed = list(seg = seg, proportionLength = FALSE))
}

bootExcludes <- function(ex) {
  list(none = list(),
       drop = list(exclude = ex),
       trim = list(exclude = ex, excludeOption = "trim"))
}

test_that("a y with reserved metadata column names fails in replicate 1", {
  inp <- invalidInputs()
  y <- inp$y
  bad <- y
  mc <- mcols(bad)
  mc$start <- start(bad) + 10L
  mc$end <- end(bad) - 10L
  ## mcols<- does not validate the GRanges
  mcols(bad) <- mc
  expect_error(validObject(bad), "names of metadata columns")
  blockLength <- 1e5
  routes <- bootRoutes(inp$seg)
  for (rn in names(routes)) for (ex in bootExcludes(inp$ex)) {
    args <- c(list(blockLength = blockLength), routes[[rn]], ex)
    got <- bootOutcome(function() do.call(bootRanges, c(list(bad, R = 3), args)))
    expect_s3_class(got$error, "simpleError")
    expect_match(conditionMessage(got$error),
                 "names of metadata columns cannot be one of", info = rn)
    expect_identical(got$warnings, character(), info = rn)
    ## the error comes from the first shift() of replicate 1: across
    ## chromosomes after its draws, which do not read metadata columns;
    ## within chromosomes after the draws for the first chromosome
    if (rn == "within") {
      set.seed(1)
      runif(ceiling(seqlengths(y)[["chr1"]] / blockLength))
      expected <- get(".Random.seed", envir = globalenv())
    } else if (rn == "permuteWithin") {
      set.seed(1)
      sample(ceiling(seqlengths(y)[["chr1"]] / blockLength))
      expected <- get(".Random.seed", envir = globalenv())
    } else {
      expected <- bootOutcome(function()
        do.call(bootRanges, c(list(y, R = 1), args)))$seed
    }
    expectSameSeed(got$seed, expected, info = rn)
  }
})

test_that("every reserved metadata column name fails as start and end do", {
  inp <- invalidInputs()
  y <- inp$y
  reserved <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                "isCircular", "start", "end", "width", "element")
  expected <- bootOutcome(function() bootRanges(y, 1e5, R = 1))$seed
  for (nm in reserved) {
    bad <- y
    mcols(bad)[[nm]] <- seq_along(bad)
    for (ex in bootExcludes(inp$ex)[c("none", "drop")]) {
      got <- bootOutcome(function()
        do.call(bootRanges, c(list(bad, 1e5, R = 2), ex)))
      expect_match(conditionMessage(got$error),
                   "names of metadata columns cannot be one of", info = nm)
      expectSameSeed(got$seed, expected, info = nm)
    }
  }
})

test_that("a state given no random block fails in the replicate that draws it", {
  ## proportionLength = FALSE: state 2 (5 of 200 bases) gets a random block
  ## only when one of the 6 bait blocks lands in it; a replicate where none
  ## does fails moving its blocks
  y <- GRanges("chr1", IRanges(seq(1, 191, by = 10), width = 5),
               seqlengths = c(chr1 = 200))
  seg <- GRanges("chr1", IRanges(c(1, 196), c(195, 200)), state = c(1, 2))
  ex <- GRanges("chr1", IRanges(50, 60), seqlengths = c(chr1 = 200))
  msg <- "length(rearranged_blocks_start) == length(random_blocks_start) is not TRUE"
  for (exArgs in bootExcludes(ex)) {
    boot <- function(R) function()
      do.call(bootRanges, c(list(y, blockLength = 40, seg = seg,
                                 proportionLength = FALSE, R = R), exArgs))
    ## replicate 1 fails: the stream stops after its draws whatever R is
    one <- bootOutcome(boot(1), seed = 1)
    expect_identical(conditionMessage(one$error), msg)
    for (R in c(3, 10)) {
      got <- bootOutcome(boot(R), seed = 1)
      expect_identical(conditionMessage(got$error), msg)
      expectSameSeed(got$seed, one$seed)
    }
    ## replicate 1 succeeds and replicate 2 fails: the stream stops after
    ## replicate 2's draws
    expect_null(bootOutcome(boot(1), seed = 6)$error)
    two <- bootOutcome(boot(2), seed = 6)
    expect_identical(conditionMessage(two$error), msg)
    five <- bootOutcome(boot(5), seed = 6)
    expect_identical(conditionMessage(five$error), msg)
    expectSameSeed(five$seed, two$seed)
  }
})

test_that("R that is not a count fails before any draw", {
  inp <- invalidInputs()
  routes <- bootRoutes(inp$seg)
  for (rn in names(routes)) {
    for (R in list(NA, -1, NULL, c(2, 3))) {
      got <- bootOutcome(function()
        do.call(bootRanges, c(list(inp$y, 1e5, R = R), routes[[rn]])))
      ## integer(R), in replicate(R, ...)
      expect_identical(conditionMessage(got$error),
                       "invalid 'length' argument", info = rn)
      expectSameSeed(got$seed, seedAfterSet(), info = rn)
    }
    ## no replicate: there are no ranges to label
    got <- bootOutcome(function()
      do.call(bootRanges, c(list(inp$y, 1e5, R = 0), routes[[rn]])))
    expect_match(conditionMessage(got$error), "mcols", info = rn)
    expectSameSeed(got$seed, seedAfterSet(), info = rn)
    ## R = 2.5 makes 2 replicates
    twoHalf <- bootOutcome(function()
      do.call(bootRanges, c(list(inp$y, 1e5, R = 2.5), routes[[rn]])))
    two <- bootOutcome(function()
      do.call(bootRanges, c(list(inp$y, 1e5, R = 2), routes[[rn]])))
    expect_identical(twoHalf$value, two$value, info = rn)
    expectSameSeed(twoHalf$seed, two$seed, info = rn)
  }
})

test_that("blockLength below 1 fails as each route's first step fails", {
  inp <- invalidInputs()
  routes <- bootRoutes(inp$seg)
  ## the block count is Inf (0) or negative (-5): sample() and runif() fail
  ## on it before drawing, tileGenome() for the permutation (before its
  ## draw), rep() for the permutation within chromosomes
  expected <- list(
    across = "invalid 'size' argument", segmented = "invalid 'size' argument",
    segmentedFixed = "invalid 'size' argument",
    permute = "'tilewidth' must be >= 1 and <= genome size",
    within = "invalid arguments", permuteWithin = "invalid 'times' argument")
  for (rn in names(routes)) for (bl in c(0, -5)) {
    got <- bootOutcome(function()
      do.call(bootRanges, c(list(inp$y, bl, R = 2), routes[[rn]])))
    expect_identical(conditionMessage(got$error), expected[[rn]],
                     info = paste(rn, bl))
    ## a block count of Inf is NA as an integer
    coerced <- bl == 0 && rn %in% c("across", "segmented", "segmentedFixed")
    expect_identical(got$warnings,
                     if (coerced) "NAs introduced by coercion to integer range"
                     else character(), info = paste(rn, bl))
    expectSameSeed(got$seed, seedAfterSet(), info = paste(rn, bl))
  }
})

test_that("blockLength 0.5 fails tiling the genome after replicate 1's draws", {
  y <- GRanges(c("chr1", "chr2"), IRanges(c(10, 20), width = 5),
               seqlengths = c(chr1 = 2000, chr2 = 1500))
  got <- bootOutcome(function() bootRanges(y, 0.5, R = 3))
  expect_identical(conditionMessage(got$error),
                   "'tilewidth' must be >= 1 and <= genome size")
  ## the bootstrap across chromosomes draws, then tiles the genome
  set.seed(1)
  L_s <- seqlengths(y)
  n_per_chrom <- ceiling(L_s / 0.5)
  n <- sum(n_per_chrom)
  random_chr <- sample(names(L_s), n, replace = TRUE, prob = L_s)
  runif(n, 1, (n_per_chrom[random_chr] - 1) * 0.5 + 1)
  expectSameSeed(got$seed, get(".Random.seed", envir = globalenv()))
})

test_that("storeBlockLength is read as a condition after the replicates", {
  inp <- invalidInputs()
  routes <- bootRoutes(inp$seg)
  for (rn in names(routes)) for (ex in bootExcludes(inp$ex)[c("none", "trim")]) {
    args <- c(list(inp$y, 1e5, R = 2), routes[[rn]], ex)
    ref <- bootOutcome(function()
      do.call(bootRanges, c(args, storeBlockLength = FALSE)))
    expected <- list(yes = "argument is not interpretable as logical",
                     na = "missing value where TRUE/FALSE needed",
                     null = "argument is of length zero")
    values <- list(yes = "yes", na = NA, null = NULL)
    for (k in names(values)) {
      got <- bootOutcome(function()
        do.call(bootRanges, c(args, list(storeBlockLength = values[[k]]))))
      expect_identical(conditionMessage(got$error), expected[[k]],
                       info = paste(rn, k))
      expectSameSeed(got$seed, ref$seed, info = paste(rn, k))
    }
  }
})

test_that("each replicate is the one-replicate call that follows the ones before", {
  ## replicate(R, ...) evaluates its body R times in turn, so replicate r of
  ## an R-replicate call is a one-replicate call made after r - 1 others
  ranges_of <- function(br) {
    list(as.character(seqnames(br)), start(br), end(br),
         as.character(strand(br)), br$score)
  }
  expectReplicatesInTurn <- function(args, info) {
    set.seed(1)
    three <- do.call(bootRanges, c(args, R = 3))
    set.seed(1)
    for (r in 1:3) {
      one <- do.call(bootRanges, c(args, R = 1))
      expect_identical(ranges_of(three[three$iter == r]), ranges_of(one),
                       info = paste(info, r))
    }
  }
  inp <- invalidInputs()
  routes <- bootRoutes(inp$seg)
  for (rn in names(routes)) {
    expectReplicatesInTurn(c(list(inp$y, 1e5), routes[[rn]]), rn)
  }
  ## a fractional blockLength without proportional lengths: the 1235-base
  ## segment needs two blocks of 1234.5 but its tiling seq() has one start,
  ## so the target sequences outnumber the target starts by one and base
  ## pairs them within each replicate; the segments alternate states and
  ## sequences, so a target paired one off changes sequence
  sl <- c(chr1 = 50000L, chr2 = 50000L)
  y <- GRanges(rep(names(sl), each = 250),
               IRanges(rep(seq(1, 49801, by = 200), 2), width = 50),
               strand = rep(c("+", "-"), 250), seqlengths = sl)
  y$score <- seq_along(y)
  cuts <- c(1L, 1236L, seq(1236L + 2469L, 50000L, by = 2469L), 50001L)
  seg <- c(GRanges("chr1", IRanges(head(cuts, -1), end = tail(cuts, -1) - 1L)),
           GRanges("chr2", IRanges(head(cuts, -1), end = tail(cuts, -1) - 1L)))
  seqlengths(seg) <- sl
  seg$state <- rep(c(1L, 2L), length.out = length(seg))
  expectReplicatesInTurn(list(y, 1234.5, seg = seg, proportionLength = FALSE),
                         "fractional blockLength")
})
