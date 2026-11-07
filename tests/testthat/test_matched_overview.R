library(nullranges)

## A hand-built Matched object: pool rows 1, 2 and 4 are the matches, so
## rows 3, 5 and 6 are unmatched. g is a factor whose level z is absent
## from focal and matched; l is logical and h character.
hand_matched <- function() {
  pool <- data.frame(a = c(2, 4, 6, 8, 10, 12),
                     l = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
                     g = factor(c("x", "y", "x", "y", "z", "z"),
                                levels = c("x", "y", "z")),
                     h = c("u", "v", "v", "u", "u", "v"),
                     ps = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.1))
  focal <- data.frame(a = c(1, 3, 5), l = c(TRUE, FALSE, TRUE),
                      g = factor(c("x", "x", "y"), levels = c("x", "y", "z")),
                      h = c("u", "v", "u"), ps = c(0.6, 0.5, 0.4))
  idx <- c(1L, 2L, 4L)
  sets <- list(focal = focal, matched = pool[idx, ], pool = pool,
               unmatched = pool[-idx, ])
  md <- do.call(rbind, lapply(names(sets), function(s) {
    cbind(id = factor(as.integer(s == "focal"), levels = c(0, 1)),
          sets[[s]], set = s)
  }))
  Matched(matchedData = data.table::as.data.table(md), matchedIndex = idx,
          covar = c("a", "l", "g", "h"), method = "rejection",
          replace = FALSE)
}

## Reference for overview(): N, then mean and sd of numeric columns or the
## count of each level of categorical columns, for each set in the order the
## sets appear; diffData is focal minus matched before rounding.
overview_ref <- function(md, digits) {
  md <- as.data.frame(md)
  cols <- setdiff(names(md), c("id", "set"))
  summarize <- function(s, round) {
    sub <- md[md$set == s, , drop = FALSE]
    unlist(lapply(setNames(cols, cols), function(cl) {
      v <- sub[[cl]]
      if (is.numeric(v)) {
        out <- c(mean = mean(v), sd = sd(v))
        if (round) signif(out, digits) else out
      } else {
        tab <- table(as.factor(v))
        setNames(as.vector(tab), names(tab))
      }
    }))
  }
  sets <- unique(md$set)
  agg <- do.call(rbind, lapply(sets, function(s) {
    c(N = sum(md$set == s), summarize(s, TRUE))
  }))
  agg <- data.frame(set = sets, agg, check.names = FALSE)
  d <- signif(summarize("focal", FALSE) - summarize("matched", FALSE), digits)
  list(aggData = agg,
       diffData = as.data.frame(as.list(d), check.names = FALSE),
       quality = abs(d[["ps.mean"]]))
}

## Same summaries, compared as plain data.frames
expect_same_overview <- function(a, b) {
  expect_identical(as.data.frame(a$aggData), as.data.frame(b$aggData))
  expect_identical(as.data.frame(a$diffData), as.data.frame(b$diffData))
  expect_identical(a$quality, b$quality)
  expect_identical(a$MatchedClass, b$MatchedClass)
}

test_that("overview of a hand-built Matched object matches hand calculation", {

  ov <- overview(hand_matched())
  expect_s4_class(ov, "MatchedOverview")

  ## Two significant digits. a: focal (1, 3, 5) mean 3, sd 2; matched
  ## (2, 4, 8) mean 4.67, sd 3.06; pool (2, ..., 12) mean 7, sd 3.74;
  ## unmatched (6, 10, 12) mean 9.33, sd 3.06. ps: focal (.6, .5, .4) mean
  ## .5, sd .1; matched (.5, .4, .2) mean .367, sd .153; pool mean .267,
  ## sd .163; unmatched (.3, .1, .1) mean .167, sd .115.
  expected <- data.frame(
    set = c("focal", "matched", "pool", "unmatched"),
    N = c(3, 3, 6, 3),
    a.mean = c(3, 4.7, 7, 9.3), a.sd = c(2, 3.1, 3.7, 3.1),
    l.FALSE = c(1, 2, 3, 1), l.TRUE = c(2, 1, 3, 2),
    g.x = c(2, 1, 2, 1), g.y = c(1, 2, 2, 0), g.z = c(0, 0, 2, 2),
    h.u = c(2, 2, 3, 1), h.v = c(1, 1, 3, 2),
    ps.mean = c(0.5, 0.37, 0.27, 0.17), ps.sd = c(0.1, 0.15, 0.16, 0.12)
  )
  expect_equal(as.data.frame(ov$aggData), expected)

  ## focal - matched from unrounded values: a.mean 3 - 4.667, a.sd 2 - 3.055,
  ## ps.mean .5 - .3667, ps.sd .1 - .1528; counts subtract.
  expected_diff <- data.frame(
    a.mean = -1.7, a.sd = -1.1, l.FALSE = -1, l.TRUE = 1,
    g.x = 1, g.y = -1, g.z = 0, h.u = 0, h.v = 0,
    ps.mean = 0.13, ps.sd = -0.053
  )
  expect_equal(as.data.frame(ov$diffData), expected_diff)
  expect_equal(ov$quality, 0.13)
  expect_identical(as.vector(ov$MatchedClass), "Matched")
})

test_that("overview digits set the significant digits of every summary", {

  x <- hand_matched()
  md <- matchedData(x)
  for (digits in list(1, 3, 4, 3L)) {
    ov <- overview(x, digits = digits)
    ref <- overview_ref(md, digits)
    expect_equal(as.data.frame(ov$aggData), ref$aggData)
    expect_equal(as.data.frame(ov$diffData), ref$diffData)
    expect_equal(ov$quality, ref$quality)
  }
  ## a.sd of the pool is sd(c(2, 4, ..., 12)) = 3.742
  expect_equal(overview(x, digits = 3)$aggData$a.sd[3], 3.74)
  expect_equal(overview(x, digits = 1)$aggData$a.sd[3], 4)

  ## digits defaults to 2; signif() treats digits below 1 as 1
  expect_same_overview(overview(x), overview(x, digits = 2))
  expect_same_overview(overview(x, digits = 0), overview(x, digits = 1))
  expect_same_overview(overview(x, digits = -1), overview(x, digits = 1))
})

test_that("overview agrees with base R on the example data set", {

  set.seed(123)
  x <- makeExampleMatchedDataSet(type = "DataFrame", matched = TRUE,
                                 method = "nearest", replace = TRUE)
  md <- matchedData(x)
  for (digits in c(2, 5)) {
    ov <- overview(x, digits = digits)
    ref <- overview_ref(md, digits)
    expect_equal(as.data.frame(ov$aggData), ref$aggData)
    expect_equal(as.data.frame(ov$diffData), ref$diffData)
    expect_equal(ov$quality, ref$quality)
  }
  ov <- overview(x)
  expect_identical(ov$aggData$set, c("focal", "matched", "pool", "unmatched"))
  expect_identical(ov$aggData$N,
                   c(500, 500, 1e4, 1e4 - length(unique(indices(x)))))
  expect_identical(names(ov$aggData),
                   c("set", "N", "feature2.mean", "feature2.sd",
                     paste0("feature3.", letters[1:5]), "ps.mean", "ps.sd"))
})

test_that("overview records the class of the Matched object", {

  p <- data.frame(a = seq(0.5, 15, by = 0.5), b = seq_len(30))
  f <- data.frame(a = c(1, 1, 2, 2, 3, 6), b = 1:6)
  gr <- GRanges("chr1", IRanges(start = 10 * seq_len(30), width = 5))
  mcols(gr) <- DataFrame(p)
  fgr <- GRanges("chr2", IRanges(start = 10 * seq_len(6), width = 5))
  mcols(fgr) <- DataFrame(f)
  reg <- GRanges("chr1", IRanges(start = seq(1, 3601, by = 100), width = 50))
  gi <- GInteractions(anchor1 = 1:30, anchor2 = 7:36, regions = reg)
  mcols(gi) <- DataFrame(p)
  fgi <- GInteractions(anchor1 = 1:6, anchor2 = 2:7, regions = reg)
  mcols(fgi) <- DataFrame(f)

  objs <- list(
    MatchedDataFrame = matchRanges(f, p, ~a, method = "nearest"),
    MatchedGRanges = matchRanges(fgr, gr, ~a, method = "nearest"),
    MatchedGInteractions = matchRanges(fgi, gi, ~a, method = "nearest")
  )
  for (cl in names(objs)) {
    ov <- overview(objs[[cl]])
    expect_identical(as.vector(ov$MatchedClass), cl)
    ref <- overview_ref(matchedData(objs[[cl]]), 2)
    expect_equal(as.data.frame(ov$aggData), ref$aggData)
    expect_equal(as.data.frame(ov$diffData), ref$diffData)
    expect_equal(ov$quality, ref$quality)
  }
  ## The three objects hold the same matchedData, so the same summaries
  ov <- lapply(objs, overview)
  for (i in 2:3) {
    expect_identical(as.data.frame(ov[[i]]$aggData),
                     as.data.frame(ov[[1]]$aggData))
    expect_identical(as.data.frame(ov[[i]]$diffData),
                     as.data.frame(ov[[1]]$diffData))
  }
})

test_that("overview with a one-row focal set and no unmatched rows", {

  ## A single focal row has no sd
  f <- data.frame(a = 3, b = 1L)
  p <- data.frame(a = c(1, 2, 3, 4, 5, 6), b = 1:6)
  ov <- overview(matchRanges(f, p, ~a, method = "nearest"))
  expect_identical(ov$aggData$N, c(1, 1, 6, 5))
  expect_equal(ov$aggData$a.mean, c(3, 3, 3.5, 3.6))
  ## sd(1:6) = 1.87; sd(c(1, 2, 4, 5, 6)) = sqrt(17.2 / 4) = 2.07
  expect_true(all(is.na(ov$aggData$a.sd[1:2])))
  expect_equal(ov$aggData$a.sd[3:4], c(1.9, 2.1))
  expect_true(all(is.na(ov$aggData$ps.sd[1:2])))
  expect_true(is.na(ov$diffData$a.sd))
  expect_identical(ov$diffData$a.mean, 0)

  ## Every pool row matched: no unmatched row in the summary
  f <- data.frame(a = c(1, 2, 3, 3, 2), b = 1:5)
  p <- data.frame(a = c(1, 2, 3), b = 1:3)
  x <- matchRanges(f, p, ~a, method = "nearest")
  ov <- overview(x)
  expect_identical(ov$aggData$set, c("focal", "matched", "pool"))
  expect_identical(ov$aggData$N, c(5, 5, 3))
  ref <- overview_ref(matchedData(x), 2)
  expect_equal(as.data.frame(ov$aggData), ref$aggData)
})

test_that("MatchedOverview $, [[ and names", {

  ov <- overview(hand_matched())
  expect_identical(names(ov),
                   c("MatchedClass", "aggData", "diffData", "quality"))
  for (s in names(ov)) {
    expect_identical(do.call(`$`, list(ov, s)), slot(ov, s))
  }
  expect_identical(ov[["quality"]], ov@quality)
  expect_identical(ov[[4]], ov@quality)
  expect_identical(as.vector(ov[["MatchedClass"]]), "Matched")
  expect_identical(as.vector(ov[[1]]), "Matched")
  expect_error(ov$notASlot, "no slot of name")
})

test_that("show prints both summary tables of a MatchedOverview", {

  ov <- overview(hand_matched())
  out <- capture.output(show(ov))
  expect_identical(out[1], "Matched object: ")
  sep <- which(out == "--------")
  expect_length(sep, 1)
  expect_identical(out[sep + 1], "focal - matched: ")
  expect_identical(out[2:(sep - 1)],
                   capture.output(print(ov$aggData, row.names = FALSE)))
  expect_identical(out[(sep + 2):length(out)],
                   capture.output(print(ov$diffData, row.names = FALSE)))

  ## Auto-printing goes through show()
  expect_identical(capture.output(print(ov)), out)

  x <- matchRanges(data.frame(a = c(1, 2), b = 1:2),
                   data.frame(a = c(1, 2, 3), b = 1:3), ~a)
  expect_identical(capture.output(show(overview(x)))[1],
                   "MatchedDataFrame object: ")
})
