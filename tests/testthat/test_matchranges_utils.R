library(nullranges)

test_that("combnCov lists every non-empty combination of covariates", {
  ## the combinations are the rows of expand.grid over {name, absent} per
  ## covariate, first covariate varying fastest, without the all-absent row
  expect_identical(combnCov("a"), "~a")
  expect_identical(combnCov(c("a", "b")), c("~a+b", "~b", "~a"))
  expect_identical(combnCov(c("a", "b", "c")),
                   c("~a+b+c", "~b+c", "~a+c", "~c", "~a+b", "~b", "~a"))

  f <- combnCov(c("x1", "x2", "x3", "x4"))
  expect_length(f, 2^4 - 1)
  expect_false(anyDuplicated(f) > 0)
  for (s in f) {
    expect_true(all(labels(terms(as.formula(s))) %in% c("x1", "x2", "x3", "x4")))
  }
})

## The example data drawn the way makeExampleMatchedDataSet draws it
example_reference <- function() {
  data.frame(
    feature1 = rep(c(TRUE, FALSE), c(500, 1e4)),
    feature2 = c(abs(rnorm(500, mean = 4, sd = 2)),
                 runif(1e4, min = 0, max = 12)),
    feature3 = c(sample(letters[1:5], 500, replace = TRUE,
                        prob = c(0.1, 0.3, 0.4, 0.1, 0.05)),
                 sample(letters[1:5], 1e4, replace = TRUE,
                        prob = c(0.4, 0.3, 0.1, 0.1, 0.05))))
}

test_that("makeExampleMatchedDataSet builds each input class from the same draws", {
  set.seed(1)
  ref <- example_reference()

  set.seed(1)
  expect_identical(makeExampleMatchedDataSet(type = "data.frame"), ref)

  set.seed(1)
  x <- makeExampleMatchedDataSet(type = "data.table")
  expect_s3_class(x, "data.table")
  expect_identical(as.data.frame(x), ref)

  ## DataFrame is the default type
  set.seed(1)
  x <- makeExampleMatchedDataSet()
  expect_s4_class(x, "DFrame")
  expect_identical(x, DataFrame(ref))

  set.seed(1)
  x <- makeExampleMatchedDataSet(type = "GRanges")
  expect_s4_class(x, "GRanges")
  expect_identical(as.character(seqnames(x)), rep("chr1", 10500))
  expect_identical(start(x), seq_len(10500))
  expect_identical(width(x), rep(100L, 10500))
  expect_identical(mcols(x), DataFrame(ref))

  set.seed(1)
  x <- makeExampleMatchedDataSet(type = "GInteractions")
  expect_s4_class(x, "GInteractions")
  expect_identical(anchors(x, "first"), granges(anchors(x, "second")))
  expect_identical(start(anchors(x, "first")), seq_len(10500))
  expect_identical(width(anchors(x, "first")), rep(100L, 10500))
  expect_identical(mcols(x), DataFrame(ref))

  expect_error(makeExampleMatchedDataSet(type = "matrix"), "should be one of")
})

test_that("makeExampleMatchedDataSet(matched = TRUE) matches feature1 on the others", {
  skip_if_not_installed("ks")
  ## rejection without replacement is the default
  set.seed(123)
  m <- makeExampleMatchedDataSet(type = "GInteractions", matched = TRUE)
  set.seed(123)
  x <- makeExampleMatchedDataSet(type = "GInteractions")
  ref <- matchRanges(x[x$feature1], x[!x$feature1],
                     covar = ~feature2 + feature3,
                     method = "rejection", replace = FALSE)
  expect_s4_class(m, "MatchedGInteractions")
  expect_identical(method(m), "rejection")
  expect_false(withReplacement(m))
  expect_identical(covariates(m), c("feature2", "feature3"))
  expect_identical(indices(m), indices(ref))
  expect_identical(focal(m), x[x$feature1])
  expect_identical(pool(m), x[!x$feature1])

  set.seed(123)
  m <- makeExampleMatchedDataSet(type = "data.frame", matched = TRUE,
                                 method = "nearest", replace = TRUE)
  set.seed(123)
  x <- makeExampleMatchedDataSet(type = "data.frame")
  ref <- matchRanges(x[x$feature1, ], x[!x$feature1, ],
                     covar = ~feature2 + feature3,
                     method = "nearest", replace = TRUE)
  expect_s4_class(m, "MatchedDataFrame")
  expect_identical(method(m), "nearest")
  expect_true(withReplacement(m))
  expect_identical(indices(m), indices(ref))
  expect_length(indices(m), 500)
})
