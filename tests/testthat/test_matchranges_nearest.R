library(nullranges)

## Propensity scores of the focal and pool sets as stored in the result
scores <- function(m) {
  md <- matchedData(m)
  list(focal = md$ps[md$set == "focal"], pool = md$ps[md$set == "pool"])
}

## Brute-force nearest neighbor matching with replacement: every pool element
## at the smallest score distance is a candidate, and the focal elements that
## share a score draw among the candidates together, in increasing score order
nearest_reference <- function(fps, pps) {
  out <- integer(length(fps))
  for (u in sort(unique(fps))) {
    d <- abs(pps - u)
    near <- which(d == min(d))
    who <- which(fps == u)
    out[who] <- near[sample(length(near), length(who), replace = TRUE)]
  }
  out
}

test_that("nearest matching picks the pool element with the closest score", {
  set.seed(1)
  focal <- data.frame(x = rnorm(25, 1), y = runif(25))
  pool <- data.frame(x = rnorm(250), y = runif(250))

  set.seed(2)
  m <- matchRanges(focal, pool, covar = ~x + y,
                   method = "nearest", replace = TRUE)
  expect_s4_class(m, "MatchedDataFrame")
  expect_identical(method(m), "nearest")
  expect_true(withReplacement(m))

  ## scores are the fitted values of a logistic regression of focal vs pool
  id <- factor(rep(c(1, 0), c(25, 250)))
  fit <- glm(id ~ x + y, family = binomial("logit"),
             data = cbind(id = id, rbind(focal, pool)))
  s <- scores(m)
  expect_equal(s$focal, unname(fitted(fit))[1:25], tolerance = 1e-8)
  expect_equal(s$pool, unname(fitted(fit))[-(1:25)], tolerance = 1e-8)

  idx <- indices(m)
  expect_length(idx, 25)
  for (i in seq_along(idx)) {
    d <- abs(s$pool - s$focal[i])
    expect_true(idx[i] %in% which(d == min(d)))
  }

  ## continuous scores have no ties, so no draw changes the result
  set.seed(2)
  expect_identical(idx, nearest_reference(s$focal, s$pool))
})

test_that("nearest matching draws among pool elements with equal scores", {
  set.seed(3)
  ## pool: 40 copies each of x = 1..4 (tied scores), then distinct values
  ## from 4.5 up, so a focal value at or above 5 never reaches a tied group
  pool <- data.frame(x = c(rep(1:4, each = 40), 0.5, 4.5,
                           sort(runif(30, 4.6, 8))))
  focal <- data.frame(x = c(sample(1:4, 20, replace = TRUE),
                            runif(4, 5, 7.5), 9, 0))

  set.seed(4)
  m <- matchRanges(focal, pool, covar = ~x, method = "nearest", replace = TRUE)
  next_draw <- runif(1)
  idx <- indices(m)
  s <- scores(m)

  ## the score is monotone in x, so equal x values match exactly and values
  ## outside the pool's range match the pool's extremes
  tied <- 1:20
  expect_identical(pool$x[idx[tied]], focal$x[tied])
  expect_identical(pool$x[idx[25]], max(pool$x))
  expect_identical(pool$x[idx[26]], 0.5)

  ## draws spread over the 40 tied elements of a level
  for (lev in unique(focal$x[tied])) {
    expect_gt(length(unique(idx[tied][focal$x[tied] == lev])), 1)
  }

  set.seed(4)
  expect_identical(idx, nearest_reference(s$focal, s$pool))
  ## one draw per distinct focal score, single candidates included
  expect_identical(runif(1), next_draw)
})

test_that("nearest matching on a factor matches within each level", {
  ## focal a/b/c = 4/2/1 and pool a/b/c = 4/8/16, so the fitted scores are
  ## the per-level focal shares 4/8, 2/10 and 1/17
  focal <- data.frame(g = factor(rep(c("a", "b", "c"), c(4, 2, 1))))
  pool <- data.frame(g = factor(rep(c("c", "b", "a"), c(16, 8, 4))))

  set.seed(5)
  m <- matchRanges(focal, pool, covar = ~g, method = "nearest", replace = TRUE)
  s <- scores(m)
  share <- c(a = 4 / 8, b = 2 / 10, c = 1 / 17)
  expect_equal(s$focal, unname(share[as.character(focal$g)]), tolerance = 1e-6)
  expect_equal(s$pool, unname(share[as.character(pool$g)]), tolerance = 1e-6)

  idx <- indices(m)
  expect_identical(pool$g[idx], focal$g)

  set.seed(5)
  expect_identical(idx, nearest_reference(s$focal, s$pool))
})
