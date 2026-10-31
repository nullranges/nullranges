library(nullranges)
library(GenomicRanges)

## bootRanges() draws and moves the replicates in chunks, and moves the
## blocks of a replicate with at least as many blocks as 'y' has ranges in
## pieces, under a budget of blocks. Moving draws nothing, so any cut gives
## the value, warnings and .Random.seed of base's code, which draws and
## moves one replicate at a time (here nullranges:::bootRanges_base).

## the value (or error message), warnings and .Random.seed after f() from
## set.seed(seed)
chunkOutcome <- function(f, seed = 1L) {
  set.seed(seed)
  warns <- character()
  value <- tryCatch(withCallingHandlers(f(), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) conditionMessage(e))
  list(value = value, warnings = warns,
       seed = sprintf("%d", get(".Random.seed", envir = globalenv())))
}

## 100 sorted ranges on two sequences, 20 of them sharing a start with a
## shorter one (y's NCList then visits them before it, unlike their index
## order), a two-state segmentation and an excluded region; with
## blockLength 1e4 each replicate has 182 blocks, more than 'y' has ranges
chunkInputs <- function() {
  set.seed(5)
  y <- GRanges(rep(c("chr1", "chr2"), each = 40),
               IRanges(c(sort(sample(990000, 40)), sort(sample(790000, 40))),
                       width = sample(50:5000, 80, replace = TRUE)),
               seqlengths = c(chr1 = 1e6 + 500, chr2 = 8e5 + 500))
  long <- y[sample(80, 20)]
  end(long) <- end(long) + 3000
  y <- sort(c(y, long))
  y$score <- seq_along(y)
  seg <- GRanges(c("chr1", "chr1", "chr2"),
                 IRanges(c(1, 500001, 1), c(5e5, 1e6 + 500, 8e5 + 500)),
                 state = c(1L, 2L, 1L), seqlengths = seqlengths(y))
  ex <- GRanges("chr1", IRanges(1e5, 2e5), seqlengths = seqlengths(y))
  list(y = y, seg = seg, ex = ex)
}

## replicates per chunk: one each, or 1, 2, 3, 1, ...
eachChunk <- function(n_rep, per_rep) as.list(rep.int(1L, n_rep))
mixedChunks <- function(n_rep, per_rep) {
  s <- integer()
  k <- 1L
  while (sum(s) < n_rep) {
    s <- c(s, min(k, n_rep - sum(s)))
    k <- k %% 3L + 1L
  }
  as.list(s)
}

test_that("chunked replicates and pieces of blocks give base's outcome", {
  inp <- chunkInputs()
  calls <- list(
    across = list(),
    drop = list(exclude = inp$ex),
    trim = list(exclude = inp$ex, excludeOption = "trim"),
    seg = list(seg = inp$seg),
    seg_drop = list(seg = inp$seg, exclude = inp$ex),
    noprop_trim = list(seg = inp$seg, proportionLength = FALSE,
                       exclude = inp$ex, excludeOption = "trim")
  )
  cuts <- list(each = list(eachChunk, 13), mixed = list(mixedChunks, 7))
  for (L in c(1e4, 2e5)) {
    for (nm in names(calls)) {
      args <- c(list(y = inp$y, blockLength = L, R = 5), calls[[nm]])
      base <- chunkOutcome(function() {
        do.call(nullranges:::bootRanges_base, args)
      })
      for (cut in names(cuts)) {
        local({
          local_mocked_bindings(boot_chunks = cuts[[cut]][[1]],
                                boot_block_budget = function() cuts[[cut]][[2]])
          got <- chunkOutcome(function() do.call(bootRanges, args))
          info <- paste(L, nm, cut)
          expect_identical(got$value, base$value, info = info)
          expect_identical(got$warnings, base$warnings, info = info)
          expect_identical(got$seed, base$seed, info = info)
        })
      }
    }
  }
})

test_that("the budget cuts the replicates into chunks of whole replicates", {
  local_mocked_bindings(boot_block_budget = function() 1000)
  expect_identical(nullranges:::boot_chunks(5L, 182), list(5L))
  expect_identical(nullranges:::boot_chunks(7L, 300),
                   list(3L, 3L, 1L))
  expect_identical(nullranges:::boot_chunks(3L, 5000), list(1L, 1L, 1L))
  expect_identical(nullranges:::boot_chunks(4L, NA_real_), list(4L))
})
