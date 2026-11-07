library(nullranges)

scores <- function(m) {
  md <- matchedData(m)
  list(focal = md$ps[md$set == "focal"], pool = md$ps[md$set == "pool"])
}

## One round of rejection sampling: accept each pool score with probability
## f(p) / (M g(p)), where f and g are kernel density estimates of the focal
## and pool scores, M is the largest f/g on a grid over the pool scores, and
## focal densities at or below 1e-3 count as zero
accept_draw <- function(fps, pps) {
  f <- ks::kde(fps)
  g <- ks::kde(pps)
  grid <- seq(from = quantile(pps, 0.001), to = quantile(pps, 0.999),
              length = 1000)
  fgrid <- predict(f, x = grid)
  ggrid <- predict(g, x = grid)
  stopifnot(all(fgrid >= 0), all(ggrid >= 0))
  M <- max(fgrid / ggrid)
  stopifnot(M <= 1e3)
  fp <- predict(f, x = pps)
  prob <- pmin(ifelse(fp > 1e-3, fp, 0) / (M * predict(g, x = pps)), 1)
  rbinom(length(pps), size = 1, prob = prob)
}

## Rejection matching: when the densities cannot be estimated, retry on
## scores with Gaussian noise of increasing sd, clamped to [0, 1], until
## enough pool elements are accepted; then draw the matches from the accepted.
## A noise level whose densities are unusable is skipped (rsMatch currently
## stops at such a level; the inputs below never reach one)
rejection_reference <- function(fps, pps, replace) {
  accept <- tryCatch(accept_draw(fps, pps), error = function(e) {
    for (sd in seq(1e-06, 0.1, length.out = 100)) {
      nf <- pmax(0, pmin(fps + rnorm(length(fps), 0, sd), 1))
      np <- pmax(0, pmin(pps + rnorm(length(pps), 0, sd), 1))
      a <- tryCatch(accept_draw(nf, np), error = function(e) NA)
      if (!anyNA(a) && sum(a) >= length(fps)) break
    }
    a
  })
  sample(which(accept == 1), length(fps), replace = replace)
}

test_that("rejection matching draws from the pool elements it accepts", {
  skip_if_not_installed("ks")
  ## with this seed both density estimates stay well above zero over the grid
  ## (other seeds can dip to -1e-15 in the tails, which rejectSample refuses)
  set.seed(12)
  focal <- data.frame(x = rnorm(30, 0.3))
  pool <- data.frame(x = rnorm(600))

  for (replace in c(TRUE, FALSE)) {
    set.seed(2)
    m <- matchRanges(focal, pool, covar = ~x,
                     method = "rejection", replace = replace)
    next_draw <- runif(1)
    expect_s4_class(m, "MatchedDataFrame")
    expect_identical(method(m), "rejection")
    expect_identical(withReplacement(m), replace)

    idx <- indices(m)
    expect_length(idx, 30)
    expect_true(all(idx >= 1 & idx <= 600))
    if (!replace) expect_false(anyDuplicated(idx) > 0)

    s <- scores(m)
    set.seed(2)
    expect_identical(idx, rejection_reference(s$focal, s$pool, replace))
    expect_identical(runif(1), next_draw)
  }
})

test_that("rejection matching retries with noisy scores when the kde fails", {
  skip_if_not_installed("ks")
  ## every focal row has the same covariate value, so all focal scores are
  ## equal and the kernel density estimate of the focal scores fails
  set.seed(3)
  focal <- data.frame(x = rep(0.5, 10))
  pool <- data.frame(x = runif(400, -2, 2))

  for (replace in c(TRUE, FALSE)) {
    set.seed(4)
    m <- matchRanges(focal, pool, covar = ~x,
                     method = "rejection", replace = replace)
    next_draw <- runif(1)
    idx <- indices(m)
    expect_length(idx, 10)
    expect_true(all(idx >= 1 & idx <= 400))
    if (!replace) expect_false(anyDuplicated(idx) > 0)

    s <- scores(m)
    expect_length(unique(s$focal), 1)
    set.seed(4)
    expect_identical(idx, rejection_reference(s$focal, s$pool, replace))
    expect_identical(runif(1), next_draw)
  }
})

test_that("rejection matching stops when too few pool elements are accepted", {
  skip_if_not_installed("ks")
  ## focal scores fill [0.45, 0.55] while pool scores spread evenly over
  ## [0, 1], so only pool scores near [0.45, 0.55] (about 15) can be accepted
  fps <- seq(0.45, 0.55, length.out = 100)
  pps <- seq(0, 1, length.out = 150)
  ## the densities are usable, so the first draw is the only one
  set.seed(5)
  expect_lt(sum(accept_draw(fps, pps)), 100)
  set.seed(5)
  expect_error(nullranges:::rsMatch(fps, pps, replace = TRUE))
})

test_that("rejectSample refuses unusable density estimates", {
  skip_if_not_installed("ks")
  ## the focal scores sit in a gap of the pool scores, so the focal density
  ## exceeds the pool density there by far more than 1e3
  expect_error(
    nullranges:::rejectSample(
      fps = seq(0.4, 0.6, length.out = 50),
      pps = c(seq(0, 0.39, length.out = 1000), seq(0.61, 1, length.out = 1000))),
    "scaling factor for density of the PS for pool is > 1e3", fixed = TRUE)

  ## scores taking two values give kde estimates that dip below zero
  expect_error(
    nullranges:::rejectSample(fps = rep(c(0.1, 0.2), c(30, 20)),
                              pps = rep(c(0.1, 0.2), c(300, 200))),
    "kernel density estimates by ks::kde are negative", fixed = TRUE)

  expect_error(nullranges:::rejectSample(fps = runif(11), pps = runif(10)),
               "focal must be <= pool for method = 'rejection'.", fixed = TRUE)

  ## a usable pair returns one 0/1 draw per pool score
  fps <- seq(0.3, 0.6, length.out = 40)
  pps <- seq(0, 1, length.out = 400)
  set.seed(6)
  accept <- nullranges:::rejectSample(fps, pps)
  set.seed(6)
  expect_identical(as.numeric(accept), as.numeric(accept_draw(fps, pps)))
  expect_true(all(accept %in% c(0, 1)))
  expect_length(accept, 400)
})
