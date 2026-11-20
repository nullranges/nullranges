library(nullranges)

## rejectSample() takes the pool's kde bandwidth from poolBandwidth() only
## while kdeVersionsPinned() is TRUE (ks 1.15.3 and KernSmooth 2.23-27), and
## otherwise lets ks::kde() pick it. These tests run each call on both paths,
## the second with kdeVersionsPinned() mocked to FALSE, and require identical
## values or errors, warnings, .Random.seed afterwards and every density
## ks::kde() returned. The densities are compared too because rbinom()'s 0/1
## draws hide a slightly different bandwidth: one from 400 grid points in
## place of dpik's 401 left every value and .Random.seed unchanged.

## value or error, every warning, .Random.seed after f() from set.seed(seed),
## and the bandwidth, grid and estimate of every ks::kde() result; with
## pinned = FALSE, poolBandwidth()'s version check is mocked to fail
rejectOutcome <- function(f, seed, pinned = TRUE) {
  kdeOriginal <- ks::kde
  kdes <- list()
  local_mocked_bindings(kde = function(x, ...) {
    fhat <- kdeOriginal(x, ...)
    kdes[[length(kdes) + 1L]] <<- fhat[c("h", "eval.points", "estimate")]
    fhat
  }, .package = "ks")
  if (!pinned) local_mocked_bindings(kdeVersionsPinned = function() FALSE)
  set.seed(seed)
  warns <- list()
  value <- withCallingHandlers(
    tryCatch(f(), error = function(e)
      list(error = conditionMessage(e), call = deparse(conditionCall(e)))),
    warning = function(w) {
      warns[[length(warns) + 1L]] <<-
        c(conditionMessage(w), deparse(conditionCall(w)))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = warns,
       seed = get(".Random.seed", envir = globalenv()), kdes = kdes)
}

expectSameOnBothPaths <- function(f, seed = 1L) {
  fast <- rejectOutcome(f, seed)
  fallback <- rejectOutcome(f, seed, pinned = FALSE)
  expect_identical(fallback, fast)
  invisible(fast)
}

## rejectSample and rsMatch (replace TRUE and FALSE) on one fps/pps pair
expectRejectSameOnBothPaths <- function(fps, pps, seed = 1L) {
  ## forced here, or a runif() argument would draw after set.seed()
  force(fps)
  force(pps)
  expectSameOnBothPaths(function() rejectSample(fps, pps), seed)
  expectSameOnBothPaths(function() rsMatch(fps, pps, replace = TRUE), seed)
  expectSameOnBothPaths(function() rsMatch(fps, pps, replace = FALSE), seed)
}

skipUnlessPinned <- function() {
  skip_if_not_installed("ks")
  skip_if_not_installed("KernSmooth")
  loadNamespace("ks")
  loadNamespace("KernSmooth")
  skip_if_not(kdeVersionsPinned(),
              "ks 1.15.3 and KernSmooth 2.23-27 not loaded: no fast path")
}

test_that("the mocked version check makes poolBandwidth decline", {
  set.seed(1)
  x <- rnorm(1000)
  fallback <- function() {
    local_mocked_bindings(kdeVersionsPinned = function() FALSE)
    poolBandwidth(x)
  }
  expect_null(fallback())
  skipUnlessPinned()
  expect_true(is.double(poolBandwidth(x)))
})

test_that("rejectSample and rsMatch match on benchmark-sized pools", {
  skipUnlessPinned()
  set.seed(11)
  fps <- plogis(rnorm(10000, -2))
  pps <- plogis(rnorm(500000, -3))
  expect_false(is.null(poolBandwidth(pps)))
  ## the focal and the pool density
  res <- expectSameOnBothPaths(function() rejectSample(fps, pps), 2L)
  expect_length(res$kdes, 2L)
  expectRejectSameOnBothPaths(fps, pps, 2L)
})

test_that("rejectSample and rsMatch match on small pools", {
  skipUnlessPinned()
  set.seed(12)
  for (np in c(2L, 3L, 10L, 60L, 400L)) {
    pps <- runif(np)
    for (nf in unique(pmin(c(1L, 2L, 5L, 50L), np))) {
      expectRejectSameOnBothPaths(runif(nf), pps, np + nf)
    }
  }
})

test_that("rejectSample and rsMatch match on tied and constant scores", {
  skipUnlessPinned()
  set.seed(13)
  ## three tied values: a positive scale estimate, so the fast path runs
  pps <- sample(c(0.1, 0.2, 0.3), 3000, replace = TRUE)
  expect_false(is.null(poolBandwidth(pps)))
  expectRejectSameOnBothPaths(sample(c(0.2, 0.3), 100, replace = TRUE), pps)
  ## mostly one value: interquartile range 0, so the fast path declines
  pps <- c(rep(0.4, 2000), runif(200))
  expect_null(poolBandwidth(pps))
  expectRejectSameOnBothPaths(runif(50, 0.3, 0.5), pps)
  ## constant pool, and constant focal and pool
  expectRejectSameOnBothPaths(runif(30), rep(0.25, 500))
  expectRejectSameOnBothPaths(rep(0.25, 30), rep(0.25, 500))
})

test_that("rejectSample and rsMatch match where the fast path declines", {
  skipUnlessPinned()
  set.seed(14)
  fps <- runif(40)
  pps <- runif(1000)
  declines <- list(
    na = c(pps, NA),
    nan = c(pps, NaN),
    inf = c(pps, Inf),
    integer = sample(1000L, 1000L, replace = TRUE),
    attribute = structure(pps, extra = 1),
    one = 0.5)
  for (p in declines) {
    expect_null(poolBandwidth(p))
    expectRejectSameOnBothPaths(fps, p)
  }
  ## names and a far outlier do not stop the fast path; the outlier makes
  ## the bandwidth's binning grid too coarse, which bkfe warns about
  named <- setNames(pps, paste0("p", seq_along(pps)))
  expect_false(is.null(poolBandwidth(named)))
  expectRejectSameOnBothPaths(fps, named)
  outlier <- c(rnorm(3000), 1e6)
  expect_false(is.null(suppressWarnings(poolBandwidth(outlier))))
  res <- expectSameOnBothPaths(function() rejectSample(rnorm(100), outlier))
  expect_true(length(res$warnings) > 0L)
  expectRejectSameOnBothPaths(rnorm(100), outlier)
})

test_that("matchRanges(method = 'rejection') matches on both paths", {
  skipUnlessPinned()
  suppressPackageStartupMessages(library(data.table))
  set.seed(15)
  makeSet <- function(n, shift) {
    data.table(a = rnorm(n, shift), b = runif(n), c = rnorm(n, shift / 2),
               g = factor(sample(c("x", "y", "z"), n, replace = TRUE)))
  }
  sizes <- list(c(10000L, 500000L), c(50L, 3000L))
  for (s in sizes) {
    focal <- makeSet(s[1L], 0.5)
    pool <- makeSet(s[2L], 0)
    for (replace in c(TRUE, FALSE)) {
      expectSameOnBothPaths(function()
        matchRanges(focal, pool, ~a + b + c, method = "rejection",
                    replace = replace), s[1L])
    }
  }
  ## a categorical covariate alone: discrete propensity scores
  focal <- makeSet(100L, 0)
  pool <- makeSet(2000L, 0)
  for (replace in c(TRUE, FALSE)) {
    expectSameOnBothPaths(function()
      matchRanges(focal, pool, ~g, method = "rejection", replace = replace))
  }
})
