library(nullranges)
library(GenomicRanges)

## matchRanges() with compiledCode = TRUE (the default) and FALSE: the value,
## the warnings and .Random.seed afterwards from the same seed are the same.

## the value (or the error and its call), warnings and .Random.seed after
## f() from set.seed(seed), the index data.table leaves on matchedData included.
matchRun <- function(f, seed = 1L) {
  set.seed(seed)
  warns <- character()
  value <- tryCatch(withCallingHandlers(f(),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(error = conditionMessage(e),
                                 call = deparse(conditionCall(e))))
  list(value = value, warnings = warns,
       seed = sprintf("%d", get(".Random.seed", envir = globalenv())))
}

test_that("compiledCode = FALSE gives the compiled result", {
  for (type in c("DataFrame", "GRanges", "GInteractions")) {
    set.seed(5)
    x <- makeExampleMatchedDataSet(type = type)
    focal <- x[x$feature1, ]
    pool <- x[!x$feature1, ]
    for (method in c("nearest", "rejection", "stratified")) {
      for (replace in c(TRUE, FALSE)) {
        if (method == "nearest" && !replace) next
        info <- paste(type, method, replace)
        run <- function(...) matchRanges(focal, pool, covar = ~feature2 + feature3,
                                         method = method, replace = replace, ...)
        fast <- matchRun(run)
        slow <- matchRun(function() run(compiledCode = FALSE))
        expect_identical(fast$value, slow$value, info = info)
        expect_identical(fast$warnings, slow$warnings, info = info)
        expect_identical(fast$seed, slow$seed, info = info)
        expect_true(is(fast$value, "Matched"), info = info)
        expect_false(is.null(attr(fast$value@matchedData, "index")), info = info)
      }
    }
  }
})

test_that("a data.frame and a one-covariate call give the same result", {
  set.seed(6)
  df <- data.frame(is_focal = rep(c(TRUE, FALSE), c(30, 200)),
                   a = c(rnorm(30, 1), rnorm(200)), b = runif(230))
  focal <- df[df$is_focal, ]
  pool <- df[!df$is_focal, ]
  for (covar in list(~a, ~a + b)) {
    fast <- matchRun(function() matchRanges(focal, pool, covar))
    slow <- matchRun(function() matchRanges(focal, pool, covar,
                                            compiledCode = FALSE))
    expect_identical(fast, slow)
  }
})

test_that("a call the compiled code declines runs the R code", {
  set.seed(7)
  x <- makeExampleMatchedDataSet(type = "DataFrame")
  focal <- x[x$feature1, ]
  pool <- x[!x$feature1, ]
  ## a covariate that is not a column: the error and its call are those of
  ## the R code
  fast <- matchRun(function() matchRanges(focal, pool, ~feature9))
  slow <- matchRun(function() matchRanges(focal, pool, ~feature9,
                                          compiledCode = FALSE))
  expect_identical(fast, slow)
  expect_match(fast$value$error, "must be columns")
  expect_match(fast$value$call, "matchRanges_Core(", fixed = TRUE)
  ## a compiled path that draws and then signals a warning or an error: the
  ## generator is put back and the value is that of the R code
  slow <- matchRun(function() matchRanges(focal, pool, ~feature2,
                                          compiledCode = FALSE))
  for (signal in list(function() warning("from the compiled path"),
                      function() stop("from the compiled path"))) {
    testthat::local_mocked_bindings(
      matchRanges_MatchedDataFrame_compiled = function(...) {
        runif(5)
        signal()
      })
    fast <- matchRun(function() matchRanges(focal, pool, ~feature2))
    expect_identical(fast, slow)
  }
})

test_that("the compiled code does not decline an ordinary call", {
  set.seed(8)
  for (type in c("DataFrame", "GRanges", "GInteractions")) {
    x <- makeExampleMatchedDataSet(type = type)
    focal <- x[x$feature1, ]
    pool <- x[!x$feature1, ]
    f <- switch(type,
      DataFrame = nullranges:::matchRanges_MatchedDataFrame_compiled,
      GRanges = nullranges:::matchRanges_MatchedGRanges_compiled,
      GInteractions = nullranges:::matchRanges_MatchedGInteractions_compiled)
    out <- nullranges:::run_or_decline(f, focal, pool, ~feature2 + feature3,
                                       "nearest", TRUE)
    expect_true(is(out, "Matched"), info = type)
  }
})
