library(nullranges)

test_that("Matched classes keep their slots and parents", {

  documented <- c(matchedData = "data.table", matchedIndex = "integer",
                  covar = "character", method = "character",
                  replace = "logical")
  expect_identical(getSlots("Matched")[names(documented)], documented)

  parents <- c(MatchedDataFrame = "DFrame",
               MatchedGRanges = "GRanges",
               MatchedGInteractions = "GInteractions")
  for (cl in names(parents)) {
    expect_true(extends(cl, "Matched"))
    expect_true(extends(cl, parents[[cl]]))
    slots <- getSlots(cl)[c("focal", "pool", "delegate")]
    expected <- if (cl == "MatchedDataFrame") "DF_OR_df_OR_dt" else parents[[cl]]
    expect_identical(unname(slots), rep(expected, 3))
  }
  expect_true(extends("BootRanges", "GRanges"))
})

test_that("Matched validity requires id, ps and set columns", {

  cols <- list(id = factor(c(1, 0, 0)), ps = c(0.6, 0.4, 0.4),
               set = c("focal", "matched", "pool"))
  md <- do.call(data.table::data.table, cols)
  x <- Matched(matchedData = md, matchedIndex = 1L, covar = "a",
               method = "nearest", replace = TRUE)
  expect_s4_class(x, "Matched")
  expect_identical(matchedData(x), md)
  expect_identical(covariates(x), "a")
  expect_identical(method(x), "nearest")
  expect_true(withReplacement(x))
  expect_identical(matchedData(MatchedDataFrame(matchedData = md)), md)
  expect_identical(matchedData(MatchedGRanges(matchedData = md)), md)
  expect_identical(matchedData(MatchedGInteractions(matchedData = md)), md)

  ## Each subclass runs the Matched validity method too
  for (drop in names(cols)) {
    bad <- do.call(data.table::data.table, cols[names(cols) != drop])
    expect_error(Matched(matchedData = bad), "matchedData")
    expect_error(MatchedDataFrame(matchedData = bad), "matchedData")
    expect_error(MatchedGRanges(matchedData = bad), "matchedData")
    expect_error(MatchedGInteractions(matchedData = bad), "matchedData")
  }
})

test_that("MatchedDataFrame takes its DataFrame slots from the delegate", {

  d <- data.frame(a = c(1.5, 2.5, 3.5), b = c("x", "y", "z"),
                  row.names = c("r1", "r2", "r3"))
  delegates <- list(data.frame = d,
                    data.table = data.table::as.data.table(d),
                    DataFrame = DataFrame(d))
  for (del in delegates) {
    x <- MatchedDataFrame(delegate = del)
    expect_s4_class(x, "MatchedDataFrame")
    expect_identical(dim(x), c(3L, 2L))
    expect_identical(colnames(x), c("a", "b"))
    expect_identical(x$a, d$a)
    expect_identical(x$b, d$b)
    expect_identical(x@delegate, DataFrame(del))
    expect_identical(rownames(x), rownames(DataFrame(del)))
  }
  expect_identical(rownames(MatchedDataFrame(delegate = d)), c("r1", "r2", "r3"))
  expect_identical(rownames(MatchedDataFrame(delegate = DataFrame(d))),
                   c("r1", "r2", "r3"))

  ## Column metadata and object metadata come along
  del <- DataFrame(d)
  mcols(del) <- DataFrame(info = c("first", "second"))
  metadata(del) <- list(source = "test")
  x <- MatchedDataFrame(delegate = del)
  expect_identical(mcols(x), mcols(del))
  expect_identical(metadata(x), list(source = "test"))
})

test_that("Matched slots pass through the subclass generators", {

  f <- DataFrame(a = c(1, 2))
  p <- DataFrame(a = c(1, 2, 3), b = c(10L, 20L, 30L))
  md <- data.table::data.table(
    id = factor(c(1, 1, 0, 0, 0, 0, 0, 0)),
    a = c(1, 2, 3, 2, 1, 2, 3, 1),
    ps = c(0.6, 0.5, 0.4, 0.5, 0.6, 0.5, 0.4, 0.6),
    set = rep(c("focal", "matched", "pool", "unmatched"), c(2, 2, 3, 1))
  )
  x <- MatchedDataFrame(focal = f, pool = p, matchedData = md,
                        matchedIndex = c(3L, 2L), covar = "a",
                        method = "rejection", replace = FALSE,
                        delegate = p[c(3, 2), ])
  expect_identical(focal(x), f)
  expect_identical(pool(x), p)
  expect_identical(matched(x), p[c(3, 2), ])
  expect_identical(unmatched(x), p[1, ])
  expect_identical(as.data.frame(x), as.data.frame(p[c(3, 2), ]))
  expect_identical(indices(x), c(3L, 2L))
  expect_identical(indices(x, "unmatched"), 1L)
  expect_identical(covariates(x), "a")
  expect_identical(method(x), "rejection")
  expect_false(withReplacement(x))
})

test_that("MatchedGRanges takes its GRanges slots from the delegate", {

  gr <- GRanges(c("chr1", "chr2", "chr1"),
                IRanges(c(10, 20, 30), width = c(5, 6, 7)),
                strand = c("+", "-", "*"),
                seqlengths = c(chr1 = 100, chr2 = 200))
  gr$score <- c(1, 2, 3)
  names(gr) <- c("g1", "g2", "g3")
  metadata(gr) <- list(source = "test")

  x <- MatchedGRanges(delegate = gr)
  expect_s4_class(x, "MatchedGRanges")
  expect_identical(x@delegate, gr)
  expect_identical(length(x), 3L)
  expect_identical(seqnames(x), seqnames(gr))
  expect_identical(ranges(x), ranges(gr))
  expect_identical(strand(x), strand(gr))
  expect_identical(seqlengths(x), c(chr1 = 100L, chr2 = 200L))
  expect_identical(names(x), c("g1", "g2", "g3"))
  expect_identical(mcols(x), mcols(gr))
  expect_identical(metadata(x), list(source = "test"))
})

test_that("MatchedGInteractions takes its slots from the delegate", {

  reg <- GRanges("chr1", IRanges(c(1, 101, 201, 301), width = 50))
  gi <- GInteractions(anchor1 = c(1, 2, 1), anchor2 = c(3, 4, 4),
                      regions = reg)
  gi$score <- c(5, 6, 7)
  names(gi) <- c("i1", "i2", "i3")
  metadata(gi) <- list(source = "test")

  x <- MatchedGInteractions(delegate = gi)
  expect_s4_class(x, "MatchedGInteractions")
  expect_identical(x@delegate, gi)
  expect_identical(length(x), 3L)
  expect_identical(anchorIds(x), anchorIds(gi))
  expect_identical(lapply(anchorIds(x), unname),
                   list(first = c(1L, 2L, 1L), second = c(3L, 4L, 4L)))
  expect_identical(regions(x), regions(gi))
  expect_identical(names(x), c("i1", "i2", "i3"))
  expect_identical(mcols(x), mcols(gi))
  expect_identical(metadata(x), list(source = "test"))
})

test_that("default delegates give empty objects", {
  x <- MatchedDataFrame()
  expect_s4_class(x, "MatchedDataFrame")
  expect_identical(dim(x), c(0L, 0L))
  expect_identical(x@delegate, DataFrame())
  x <- MatchedGRanges()
  expect_s4_class(x, "MatchedGRanges")
  expect_identical(length(x), 0L)
  expect_identical(x@delegate, GRanges())
  x <- MatchedGInteractions()
  expect_s4_class(x, "MatchedGInteractions")
  expect_identical(length(x), 0L)
  expect_identical(x@delegate, GInteractions())
})

test_that("subclasses print as their parent class with the matched set", {

  p <- data.frame(a = seq(0.5, 15, by = 0.5), b = seq_len(30))
  f <- data.frame(a = c(1, 1, 2, 2, 3, 6))
  mdf <- matchRanges(f, p, ~a, method = "nearest", replace = TRUE)
  out <- capture.output(show(mdf))
  expect_match(out[1], "^MatchedDataFrame with 6 rows and 2 columns")

  gr <- GRanges("chr1", IRanges(start = 10 * seq_len(30), width = 5))
  gr$a <- p$a
  fgr <- GRanges("chr2", IRanges(start = 10 * seq_len(6), width = 5))
  fgr$a <- f$a
  mgr <- matchRanges(fgr, gr, ~a, method = "nearest", replace = TRUE)
  out <- capture.output(show(mgr))
  expect_match(out[1], "^MatchedGRanges object with 6 ranges and 1 metadata column")

  reg <- GRanges("chr1", IRanges(start = seq(1, 3601, by = 100), width = 50))
  gi <- GInteractions(anchor1 = 1:30, anchor2 = 7:36, regions = reg)
  gi$a <- p$a
  fgi <- GInteractions(anchor1 = 1:6, anchor2 = 2:7, regions = reg)
  fgi$a <- f$a
  mgi <- matchRanges(fgi, gi, ~a, method = "nearest", replace = TRUE)
  out <- capture.output(show(mgi))
  expect_match(out[1],
               "^MatchedGInteractions object with 6 interactions and 1 metadata column")
})

test_that("BootRanges is a GRanges subclass that bootRanges returns", {

  gr <- GRanges(c("chr1", "chr1", "chr2"), IRanges(c(1, 51, 11), width = 10),
                seqlengths = c(chr1 = 100, chr2 = 100))
  gr$iter <- factor(c(1, 1, 2))
  b <- BootRanges(gr)
  expect_s4_class(b, "BootRanges")
  expect_s4_class(b, "GRanges")
  expect_identical(as(b, "GRanges"), gr)

  set.seed(1)
  bb <- bootRanges(granges(gr), blockLength = 20, R = 2)
  expect_s4_class(bb, "BootRanges")
})
