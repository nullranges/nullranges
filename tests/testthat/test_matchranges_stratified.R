library(nullranges)

scores <- function(m) {
  md <- matchedData(m)
  list(focal = md$ps[md$set == "focal"], pool = md$ps[md$set == "pool"])
}

## Iterative stratified matching. Each pass cuts the range of the unmatched
## scores into n equal-width bins (n = number of distinct unmatched scores),
## halving n until some bin holds at least one focal score and no more focal
## than pool scores. Every such bin, in order of first appearance among the
## pool scores, draws its pool matches; matched pool elements leave the pool.
## Once n has been halved it keeps halving on each later pass.
stratified_reference <- function(fps, pps, replace) {
  fo <- seq_along(fps)
  po <- seq_along(pps)
  rf <- integer()
  rp <- integer()
  halving <- FALSE
  while (length(rf) < length(fps)) {
    n <- if (halving) floor(n / 2) else length(unique(c(fps[fo], pps[po])))
    repeat {
      v <- c(fps[fo], pps[po])
      if (n >= 1 && max(v) > min(v)) {
        br <- seq(min(v), max(v), by = (max(v) - min(v)) / n)
        br[length(br)] <- max(v)
        fb <- .bincode(fps[fo], br, include.lowest = TRUE)
        pb <- .bincode(pps[po], br, include.lowest = TRUE)
      } else {
        fb <- rep(1L, length(fo))
        pb <- rep(1L, length(po))
      }
      bins <- unique(pb)
      ok <- vapply(bins, function(b) {
        k <- sum(fb == b)
        k > 0 && k <= sum(pb == b)
      }, logical(1))
      if (any(ok)) break
      halving <- TRUE
      n <- floor(n / 2)
    }
    for (b in bins[ok]) {
      f <- fo[fb == b]
      p <- po[pb == b]
      rf <- c(rf, f)
      rp <- c(rp, p[sample(length(p), length(f), replace = replace)])
    }
    fo <- setdiff(fo, rf)
    po <- setdiff(po, rp)
  }
  rp[order(rf)]
}

test_that("stratified matching halves its bins until every focal is matched", {
  ## focal a/b = 2/3 and pool a/b = 1/5 give scores 2/3 (a) and 3/8 (b).
  ## Pass 1 has two bins: b (3 focal, 5 pool) is matched, a (2 focal, 1 pool)
  ## is not. Pass 2 still finds 2 focal and 1 pool in the a bin, so it halves
  ## to one bin and the two a focal draw from all pool elements left.
  ## The scores are within a factor of 2, so the bin edges are exact.
  focal <- data.frame(g = c("a", "b", "b", "a", "b"))
  pool <- data.frame(g = c("b", "a", "b", "b", "b", "b"))

  for (replace in c(TRUE, FALSE)) {
    for (seed in 1:3) {
      set.seed(seed)
      m <- suppressMessages(matchRanges(focal, pool, covar = ~g,
                                        method = "stratified",
                                        replace = replace))
      expect_s4_class(m, "MatchedDataFrame")
      expect_identical(method(m), "stratified")
      expect_identical(withReplacement(m), replace)

      s <- scores(m)
      expect_equal(s$focal, c(2 / 3, 3 / 8, 3 / 8, 2 / 3, 3 / 8),
                   tolerance = 1e-6)

      set.seed(seed)
      b <- c(1, 3, 4, 5, 6)[sample(5, 3, replace = replace)]
      rest <- setdiff(1:6, b)
      a <- rest[sample(length(rest), 2, replace = replace)]
      expect_identical(indices(m), as.integer(c(a[1], b[1:2], a[2], b[3])))
    }
  }
})

test_that("stratified matching agrees with a reference on continuous scores", {
  set.seed(1)
  focal <- data.frame(x = rnorm(40, 0.5), y = runif(40))
  pool <- data.frame(x = rnorm(400), y = runif(400))

  for (replace in c(TRUE, FALSE)) {
    set.seed(2)
    m <- suppressMessages(matchRanges(focal, pool, covar = ~x + y,
                                      method = "stratified", replace = replace))
    next_draw <- runif(1)
    idx <- indices(m)
    expect_length(idx, 40)
    expect_true(all(idx >= 1 & idx <= 400))
    if (!replace) expect_false(anyDuplicated(idx) > 0)

    s <- scores(m)
    set.seed(2)
    expect_identical(idx, stratified_reference(s$focal, s$pool, replace))
    expect_identical(runif(1), next_draw)
  }
})

test_that("stratified matching keeps factor levels apart when bins separate them", {
  ## focal a/b/c = 3/2/1 and pool a/b/c = 3/6/9 give scores 1/2, 1/4 and 1/10;
  ## with 3 bins over [0.1, 0.5] each level falls in its own bin, and each bin
  ## has at least as many pool as focal elements, so one pass matches all
  focal <- data.frame(g = rep(c("a", "b", "c"), c(3, 2, 1)))
  pool <- data.frame(g = rep(c("c", "b", "a"), c(9, 6, 3)))

  set.seed(3)
  m <- suppressMessages(matchRanges(focal, pool, covar = ~g,
                                    method = "stratified", replace = FALSE))
  idx <- indices(m)
  expect_identical(pool$g[idx], focal$g)
  expect_false(anyDuplicated(idx) > 0)

  s <- scores(m)
  expect_equal(s$focal, c(1 / 2, 1 / 2, 1 / 2, 1 / 4, 1 / 4, 1 / 10),
               tolerance = 1e-6)
  set.seed(3)
  expect_identical(idx, stratified_reference(s$focal, s$pool, FALSE))
})
