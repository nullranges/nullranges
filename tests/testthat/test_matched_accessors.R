library(nullranges)

## Rows (or ranges) i of a DataFrame, GRanges or GInteractions
take_rows <- function(obj, i) {
  if (is(obj, "DataFrame")) obj[i, , drop = FALSE] else obj[i]
}

## Laws every Matched object from matchRanges() obeys, whatever the method:
## f and p are the focal and pool given to matchRanges()
expect_matched_laws <- function(x, f, p, covar, meth, replace) {
  if (is.data.frame(f)) {
    f <- DataFrame(f)
    p <- DataFrame(p)
  }
  nf <- NROW(f)
  np <- NROW(p)
  idx <- indices(x)
  unm <- setdiff(seq_len(np), idx)

  expect_type(idx, "integer")
  expect_length(idx, nf)
  expect_true(all(idx >= 1L & idx <= np))
  if (!replace) expect_identical(anyDuplicated(idx), 0L)
  expect_identical(indices(x, "focal"), seq_len(nf))
  expect_identical(indices(x, "pool"), seq_len(np))
  expect_identical(indices(x, "unmatched"), unm)
  expect_identical(covariates(x), covar)
  expect_identical(method(x), meth)
  expect_identical(withReplacement(x), replace)

  ## Sets as objects of the input class
  expect_identical(focal(x), f)
  expect_identical(pool(x), p)
  ## Contents only: matched() of a data.frame input keeps its row names
  if (is(p, "DataFrame")) {
    expect_identical(as.list(matched(x)), as.list(take_rows(p, idx)))
  } else {
    expect_identical(matched(x), take_rows(p, idx))
  }
  expect_identical(NROW(matched(x)), nf)
  expect_identical(unmatched(x), take_rows(p, unm))
  expect_identical(NROW(x), nf)

  ## matchedData: focal, matched, pool, unmatched rows in that order
  md <- as.data.frame(matchedData(x))
  expect_identical(names(md), c("id", covar, "ps", "set"))
  expect_identical(md$set, rep(c("focal", "matched", "pool", "unmatched"),
                               c(nf, nf, np, length(unm))))
  expect_identical(as.character(md$id),
                   rep(c("1", "0"), c(nf, nf + np + length(unm))))
  rows <- function(s) {
    r <- md[md$set == s, c(covar, "ps"), drop = FALSE]
    rownames(r) <- NULL
    r
  }
  pool_rows <- rows("pool")
  at <- function(i) {
    r <- pool_rows[i, , drop = FALSE]
    rownames(r) <- NULL
    r
  }
  expect_identical(rows("matched"), at(idx))
  expect_identical(rows("unmatched"), at(unm))
  cols <- if (is(f, "DataFrame")) f else mcols(f)
  pcols <- if (is(p, "DataFrame")) p else mcols(p)
  for (cv in covar) {
    expect_identical(rows("focal")[[cv]], as.vector(cols[[cv]]))
    expect_identical(pool_rows[[cv]], as.vector(pcols[[cv]]))
  }
}

## Pool row k has a = k / 2 and b = k. Each focal value of a occurs once in
## the pool, so nearest matching on ~a picks the pool row with the same a:
## a = 1 -> row 2, a = 2 -> row 4, a = 3 -> row 6, a = 6 -> row 12.
pool_df <- data.frame(a = seq(0.5, 15, by = 0.5), b = seq_len(30))
focal_df <- data.frame(a = c(1, 1, 2, 2, 3, 6), b = 101:106)
hand_idx <- c(2L, 2L, 4L, 4L, 6L, 12L)
hand_unm <- setdiff(seq_len(30), hand_idx)

test_that("hand-derived nearest matches for data.frame, data.table and DataFrame", {

  inputs <- list(
    data.frame = list(f = focal_df, p = pool_df),
    data.table = list(f = data.table::as.data.table(focal_df),
                      p = data.table::as.data.table(pool_df)),
    DataFrame = list(f = DataFrame(focal_df), p = DataFrame(pool_df))
  )
  for (nm in names(inputs)) {
    f <- inputs[[nm]]$f
    p <- inputs[[nm]]$p
    x <- matchRanges(f, p, ~a, method = "nearest", replace = TRUE)
    expect_s4_class(x, "MatchedDataFrame")
    expect_identical(indices(x), hand_idx)
    expect_identical(indices(x, "unmatched"), hand_unm)
    expect_identical(matched(x)$b, hand_idx)
    expect_identical(unmatched(x)$b, hand_unm)
    expect_identical(x$a, c(1, 1, 2, 2, 3, 6))
    expect_identical(focal(x)$b, 101:106)
    expect_matched_laws(x, f, p, "a", "nearest", TRUE)
  }
  ## Inputs without row names give a matched set identical to the pool rows
  x <- matchRanges(inputs$DataFrame$f, inputs$DataFrame$p, ~a,
                   method = "nearest", replace = TRUE)
  expect_identical(matched(x), pool(x)[hand_idx, ])
})

test_that("hand-derived nearest matches for GRanges", {

  p <- GRanges("chr1", IRanges(start = 10 * seq_len(30), width = 5))
  mcols(p) <- DataFrame(pool_df)
  f <- GRanges("chr2", IRanges(start = 10 * seq_len(6), width = 5))
  mcols(f) <- DataFrame(focal_df)

  x <- matchRanges(f, p, ~a, method = "nearest", replace = TRUE)
  expect_s4_class(x, "MatchedGRanges")
  expect_identical(indices(x), hand_idx)
  expect_identical(matched(x), p[hand_idx])
  expect_identical(unmatched(x), p[hand_unm])
  expect_identical(focal(x), f)
  expect_identical(pool(x), p)
  expect_identical(start(x), 10L * hand_idx)
  expect_identical(as.character(seqnames(x)), rep("chr1", 6))
  expect_identical(mcols(x), mcols(p)[hand_idx, ])
  expect_matched_laws(x, f, p, "a", "nearest", TRUE)
})

test_that("hand-derived nearest matches for GInteractions", {

  reg <- GRanges("chr1", IRanges(start = seq(1, 3601, by = 100), width = 50))
  p <- GInteractions(anchor1 = 1:30, anchor2 = 7:36, regions = reg)
  mcols(p) <- DataFrame(pool_df)
  f <- GInteractions(anchor1 = 1:6, anchor2 = 2:7, regions = reg)
  mcols(f) <- DataFrame(focal_df)

  x <- matchRanges(f, p, ~a, method = "nearest", replace = TRUE)
  expect_s4_class(x, "MatchedGInteractions")
  expect_identical(indices(x), hand_idx)
  expect_identical(matched(x), p[hand_idx])
  expect_identical(unmatched(x), p[hand_unm])
  expect_identical(focal(x), f)
  expect_identical(pool(x), p)
  expect_identical(unname(anchorIds(x)$first), hand_idx)
  expect_identical(unname(anchorIds(x)$second), hand_idx + 6L)
  expect_identical(mcols(x), mcols(p)[hand_idx, ])
  expect_matched_laws(x, f, p, "a", "nearest", TRUE)
})

test_that("matchedData holds the inputs' covariates and their propensity scores", {

  x <- matchRanges(focal_df, pool_df, ~a, method = "nearest", replace = TRUE)
  md <- as.data.frame(matchedData(x))
  expect_s3_class(matchedData(x), "data.table")
  expect_identical(levels(md$id), c("0", "1"))
  expect_identical(md$a, c(focal_df$a, pool_df$a[hand_idx], pool_df$a,
                           pool_df$a[hand_unm]))

  ## Propensity scores are the fitted values of a logistic regression of
  ## focal (1) versus pool (0) membership on the covariates
  fit <- glm(y ~ a, family = binomial,
             data = data.frame(y = rep(1:0, c(6, 30)),
                               a = c(focal_df$a, pool_df$a)))
  ps <- unname(fitted(fit))
  expect_equal(md$ps[md$set == "focal"], ps[1:6], tolerance = 1e-8)
  expect_equal(md$ps[md$set == "pool"], ps[-(1:6)], tolerance = 1e-8)
  expect_identical(md$ps[md$set == "matched"],
                   md$ps[md$set == "pool"][hand_idx])
  expect_identical(md$ps[md$set == "unmatched"],
                   md$ps[md$set == "pool"][hand_unm])
})

test_that("accessors report the arguments matchRanges was given", {

  f <- data.frame(a = c(1, 2, 3), cc = c(2, 1, 3))
  p <- data.frame(a = seq(0.5, 15, by = 0.5), cc = rep(c(1, 2, 3), 10))

  ## Generic defaults: method = 'nearest', replace = TRUE
  x <- matchRanges(f, p, ~a + cc)
  expect_identical(method(x), "nearest")
  expect_true(withReplacement(x))
  expect_identical(covariates(x), c("a", "cc"))
  expect_matched_laws(x, f, p, c("a", "cc"), "nearest", TRUE)

  ## A partial method name is stored in full; formula order is kept
  x <- matchRanges(f, p, ~cc + a, method = "n", replace = TRUE)
  expect_identical(method(x), "nearest")
  expect_identical(covariates(x), c("cc", "a"))
  expect_identical(names(matchedData(x)), c("id", "cc", "a", "ps", "set"))
})

test_that("every method and replace value gives sets that partition the pool", {

  set.seed(11)
  f <- data.frame(a = rnorm(60, mean = 3, sd = 1),
                  g = sample(c("u", "v"), 60, replace = TRUE))
  p <- data.frame(a = rnorm(1200, mean = 4, sd = 2),
                  g = sample(c("u", "v"), 1200, replace = TRUE))

  runs <- list(c("nearest", TRUE), c("stratified", TRUE),
               c("stratified", FALSE), c("rejection", TRUE),
               c("rejection", FALSE))
  for (run in runs) {
    meth <- run[1]
    replace <- as.logical(run[2])
    if (meth == "rejection") skip_if_not_installed("ks")
    set.seed(2)
    capture.output(
      x <- matchRanges(f, p, ~a + g, method = meth, replace = replace),
      type = "message"
    )
    expect_matched_laws(x, f, p, c("a", "g"), meth, replace)
  }
})

test_that("GRanges and GInteractions sets partition the pool for every method", {

  set.seed(12)
  n_p <- 400
  reg <- GRanges("chr1", IRanges(start = seq(1, by = 100, length.out = n_p + 10),
                                 width = 50))
  p_gr <- GRanges("chr1", IRanges(start = 100 * seq_len(n_p), width = 20))
  p_gr$a <- runif(n_p, 0, 10)
  f_gr <- GRanges("chr2", IRanges(start = 100 * seq_len(40), width = 20))
  f_gr$a <- rnorm(40, mean = 4, sd = 1)
  p_gi <- GInteractions(anchor1 = seq_len(n_p), anchor2 = seq_len(n_p) + 10L,
                        regions = reg)
  p_gi$a <- p_gr$a
  f_gi <- GInteractions(anchor1 = 1:40, anchor2 = 3:42, regions = reg)
  f_gi$a <- f_gr$a

  runs <- list(c("nearest", TRUE), c("stratified", FALSE),
               c("rejection", FALSE))
  for (run in runs) {
    meth <- run[1]
    replace <- as.logical(run[2])
    if (meth == "rejection") skip_if_not_installed("ks")
    set.seed(3)
    capture.output(
      xg <- matchRanges(f_gr, p_gr, ~a, method = meth, replace = replace),
      type = "message"
    )
    expect_s4_class(xg, "MatchedGRanges")
    expect_matched_laws(xg, f_gr, p_gr, "a", meth, replace)
    expect_identical(granges(xg), granges(p_gr[indices(xg)]))

    set.seed(3)
    capture.output(
      xi <- matchRanges(f_gi, p_gi, ~a, method = meth, replace = replace),
      type = "message"
    )
    expect_s4_class(xi, "MatchedGInteractions")
    expect_matched_laws(xi, f_gi, p_gi, "a", meth, replace)
    expect_identical(anchorIds(xi), anchorIds(p_gi[indices(xi)]))
  }
})

test_that("an empty unmatched set when every pool row is matched", {

  ## Each focal a equals one pool a, and every pool row is someone's match
  f <- data.frame(a = c(1, 2, 3, 3, 2), b = 1:5)
  p <- data.frame(a = c(1, 2, 3), b = c(10L, 20L, 30L))
  x <- matchRanges(f, p, ~a, method = "nearest", replace = TRUE)
  expect_identical(indices(x), c(1L, 2L, 3L, 3L, 2L))
  expect_identical(indices(x, "unmatched"), integer(0))
  expect_identical(dim(unmatched(x)), c(0L, 2L))
  expect_identical(matched(x)$b, c(10L, 20L, 30L, 30L, 20L))
  expect_false("unmatched" %in% matchedData(x)$set)
  expect_matched_laws(x, f, p, "a", "nearest", TRUE)
})

test_that("indices on a hand-built Matched follow the set sizes in matchedData", {

  md <- data.table::data.table(
    id = factor(c(1, 1, rep(0, 8))),
    a = c(1, 3, 8, 4, 2, 4, 6, 8, 2, 6),
    ps = c(0.6, 0.5, 0.2, 0.4, 0.5, 0.4, 0.3, 0.2, 0.5, 0.3),
    set = rep(c("focal", "matched", "pool", "unmatched"), c(2, 2, 4, 2))
  )
  x <- Matched(matchedData = md, matchedIndex = c(4L, 2L), covar = "a",
               method = "nearest", replace = TRUE)
  expect_identical(indices(x), c(4L, 2L))
  expect_identical(indices(x, "matched"), c(4L, 2L))
  expect_identical(indices(x, "focal"), 1:2)
  expect_identical(indices(x, "pool"), 1:4)
  expect_identical(indices(x, "unmatched"), c(1L, 3L))

  ## Set names may be abbreviated
  expect_identical(indices(x, "un"), c(1L, 3L))
  expect_identical(indices(x, "f"), 1:2)
  expect_identical(indices(x, "p"), 1:4)
  expect_error(indices(x, "all"), "should be one of")

  ## A pool row matched twice is not unmatched
  x@matchedIndex <- c(4L, 4L)
  expect_identical(indices(x, "unmatched"), 1:3)
})
