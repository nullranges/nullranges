library(nullranges)
library(GenomicRanges)

## bootRanges() with compiledCode = TRUE (the default) and FALSE: the value,
## the warnings and .Random.seed afterwards from the same seed are the same.

## the value (or the error and its call), warnings and .Random.seed after
## f() from set.seed(seed)
bootRun <- function(f, seed = 1L) {
  set.seed(seed)
  warns <- character()
  value <- tryCatch(withCallingHandlers(f(), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e),
                               call = deparse(conditionCall(e))))
  list(value = value, warnings = warns,
       seed = sprintf("%d", get(".Random.seed", envir = globalenv())))
}

bootInputs <- function() {
  set.seed(3)
  sl <- c(chr1 = 200000L, chr2 = 150000L)
  n <- 300
  chr <- sample(names(sl), n, replace = TRUE)
  start <- vapply(chr, function(s) sample.int(sl[[s]] - 400L, 1L), 1L)
  y <- sort(GRanges(chr, IRanges(start, width = sample(20:300, n, TRUE)),
                    seqlengths = sl, score = rnorm(n)))
  b <- seq(1L, 200000L, length.out = 5)
  seg <- GRanges(rep(names(sl), each = 4),
                 IRanges(c(b[-5], round(b[-5] * 0.75)),
                         width = c(diff(b), round(diff(b) * 0.75))),
                 seqlengths = sl, state = c(1, 2, 3, 1, 2, 1, 3, 2))
  seg <- sort(seg)
  ex <- sort(GRanges("chr1", IRanges(c(1000, 90000), width = 5000),
                     seqlengths = sl))
  list(y = y, seg = seg, ex = ex)
}

test_that("compiledCode = FALSE gives the compiled result", {
  inp <- bootInputs()
  calls <- list(
    unseg = list(),
    drop = list(exclude = inp$ex),
    trim = list(exclude = inp$ex, excludeOption = "trim"),
    permute = list(type = "permute"),
    within = list(withinChrom = TRUE),
    block_length = list(storeBlockLength = TRUE),
    seg = list(seg = inp$seg),
    seg_noprop = list(seg = inp$seg, proportionLength = FALSE),
    seg_drop = list(seg = inp$seg, exclude = inp$ex)
  )
  for (nm in names(calls)) {
    args <- c(list(y = inp$y, blockLength = 20000, R = 3), calls[[nm]])
    fast <- bootRun(function() do.call(bootRanges, args))
    slow <- bootRun(function()
      do.call(bootRanges, c(args, list(compiledCode = FALSE))))
    expect_identical(fast$value, slow$value, info = nm)
    expect_identical(fast$warnings, slow$warnings, info = nm)
    expect_identical(fast$seed, slow$seed, info = nm)
  }
})

test_that("a call the compiled code declines runs the R code", {
  inp <- bootInputs()
  ## an unsorted y: the error and its call are those of the R code
  y <- inp$y[c(2, 1, seq_along(inp$y)[-(1:2)])]
  fast <- bootRun(function() bootRanges(y, 20000, R = 2))
  slow <- bootRun(function() bootRanges(y, 20000, R = 2, compiledCode = FALSE))
  expect_identical(fast, slow)
  expect_match(fast$value$error, "is not TRUE")
  ## a compiled path that draws and then signals a warning or an error: the
  ## generator is put back and the value is that of the R code
  slow <- bootRun(function()
    bootRanges(inp$y, 20000, R = 2, compiledCode = FALSE))
  for (signal in list(function() warning("from the compiled path"),
                      function() stop("from the compiled path"))) {
    testthat::local_mocked_bindings(boot_ranges_run = function(...) {
      runif(5)
      signal()
    })
    fast <- bootRun(function() bootRanges(inp$y, 20000, R = 2))
    expect_identical(fast, slow)
  }
})

test_that("the compiled code does not decline an ordinary call", {
  inp <- bootInputs()
  expect_s4_class(nullranges:::bootRanges_compiled(
    inp$y, 20000, 2, NULL, NULL, "drop", TRUE, "bootstrap", FALSE, FALSE),
    "BootRanges")
  expect_s4_class(nullranges:::bootRanges_compiled(
    inp$y, 20000, 2, inp$seg, inp$ex, "drop", TRUE, "bootstrap", FALSE, FALSE),
    "BootRanges")
})
