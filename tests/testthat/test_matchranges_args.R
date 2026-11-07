library(nullranges)

## 6 focal and 30 pool rows with a numeric and a character covariate
args_data <- function(n_focal = 6, n_pool = 30) {
  set.seed(1)
  list(focal = data.frame(x = rnorm(n_focal, 1),
                          g = rep(c("a", "b"), length.out = n_focal)),
       pool = data.frame(x = rnorm(n_pool),
                         g = rep(c("a", "b", "c"), length.out = n_pool)))
}

test_that("covar must be a non-empty right-hand-side formula", {
  d <- args_data()

  expect_error(matchRanges(d$focal, d$pool, covar = y ~ x),
               "covar only accepts rhs of formula.", fixed = TRUE)

  ## formula(NULL) is a formula object of length 0
  expect_error(matchRanges(d$focal, d$pool, covar = formula(NULL)),
               "covar cannot be an empty formula.", fixed = TRUE)
})

test_that("every covariate must be a column of both focal and pool", {
  d <- args_data()
  msg <- "All variables in covar must be columns in both focal and pool."

  focal <- d$focal
  focal$z <- 1
  expect_error(matchRanges(focal, d$pool, covar = ~x + z), msg, fixed = TRUE)

  pool <- d$pool
  pool$z <- 1
  expect_error(matchRanges(d$focal, pool, covar = ~x + z), msg, fixed = TRUE)
  expect_error(matchRanges(d$focal, d$pool, covar = ~absent), msg, fixed = TRUE)
})

test_that("method defaults to nearest with replacement and accepts abbreviations", {
  d <- args_data()

  set.seed(2)
  m <- matchRanges(d$focal, d$pool, covar = ~x)
  expect_identical(method(m), "nearest")
  expect_true(withReplacement(m))

  set.seed(2)
  expect_identical(
    indices(matchRanges(d$focal, d$pool, covar = ~x, method = "n")),
    indices(m))

  set.seed(2)
  m <- suppressMessages(
    matchRanges(d$focal, d$pool, covar = ~x, method = "s", replace = FALSE))
  expect_identical(method(m), "stratified")
  expect_false(withReplacement(m))

  expect_error(matchRanges(d$focal, d$pool, covar = ~x, method = "optimal"),
               "should be one of")

  skip_if_not_installed("ks")
  set.seed(2)
  m <- matchRanges(d$focal, d$pool, covar = ~x, method = "r", replace = FALSE)
  expect_identical(method(m), "rejection")
  expect_false(withReplacement(m))
})

test_that("matching without replacement needs fewer focal than pool rows", {
  msg <- "focal must be <= pool when replace = FALSE."

  ## equal sizes are refused too
  for (n_pool in c(6, 5)) {
    d <- args_data(n_focal = 6, n_pool = n_pool)
    for (meth in c("nearest", "rejection", "stratified")) {
      expect_error(matchRanges(d$focal, d$pool, covar = ~x,
                               method = meth, replace = FALSE),
                   msg, fixed = TRUE)
    }
  }
})

test_that("nearest neighbor matching without replacement is refused", {
  d <- args_data()
  expect_error(matchRanges(d$focal, d$pool, covar = ~x,
                           method = "nearest", replace = FALSE),
               "nearest neighbor matching without replacement not available.",
               fixed = TRUE)
})

test_that("stratified matching needs fewer focal than pool rows with replacement", {
  d <- args_data(n_focal = 6, n_pool = 6)
  expect_error(matchRanges(d$focal, d$pool, covar = ~x,
                           method = "stratified", replace = TRUE),
               "focal must be <= pool for stratified sampling.", fixed = TRUE)
})

test_that("rejection matching with replacement needs no more focal than pool rows", {
  skip_if_not_installed("ks")
  d <- args_data(n_focal = 6, n_pool = 5)
  expect_error(matchRanges(d$focal, d$pool, covar = ~x,
                           method = "rejection", replace = TRUE),
               "focal must be <= pool for method = 'rejection'.", fixed = TRUE)
})
