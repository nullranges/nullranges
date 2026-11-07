library(nullranges)

## 10 treated and 40 control rows; default row names, which matchitToMatched
## reads as row positions
matchit_data <- function() {
  set.seed(1)
  data.frame(treat = rep(c(TRUE, FALSE), c(10, 40)),
             x = c(rnorm(10, 1), rnorm(40)),
             g = sample(c("u", "v"), 50, replace = TRUE))
}

test_that("matchitToMatched maps MatchIt's matched controls into the pool", {
  skip_if_not_installed("MatchIt")
  d <- matchit_data()
  m <- MatchIt::matchit(treat ~ x + g, data = d,
                        method = "nearest", replace = FALSE)
  mm <- matchitToMatched(m)

  expect_s4_class(mm, "MatchedDataFrame")
  expect_identical(method(mm), "MatchIt_nearest")
  expect_false(withReplacement(mm))
  expect_identical(covariates(mm), c("x", "g"))

  ## match.matrix names the control matched to each treated row; indices are
  ## those controls' positions among the control rows (rows 11:50)
  controls <- as.integer(m$match.matrix[, 1])
  idx <- indices(mm)
  expect_identical(sort(idx), sort(match(controls, which(!d$treat))))

  expect_identical(focal(mm), DataFrame(d[d$treat, ]))
  expect_identical(pool(mm), DataFrame(d[!d$treat, ]))
  expect_identical(matched(mm)$x, d$x[10L + idx])
  expect_identical(matched(mm)$g, d$g[10L + idx])

  md <- as.data.frame(matchedData(mm))
  sets <- c("focal", "matched", "pool")
  md <- md[md$set %in% sets, ]
  expect_identical(md$set, rep(sets, c(10, 10, 40)))
  rows <- c(1:10, 10L + idx, 11:50)
  expect_identical(md$x, d$x[rows])
  expect_equal(md$ps, unname(m$distance)[rows], tolerance = 1e-8)
})

test_that("matchitToMatched attaches supplied ranges with or without their mcols", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("mariner")
  d <- matchit_data()
  m <- MatchIt::matchit(treat ~ x + g, data = d,
                        method = "nearest", replace = FALSE)
  pos <- as.integer(m$match.matrix[, 1]) - 10L

  gr <- GRanges("chr2", IRanges(start = 1:50 * 100, width = 20), score = 50:1)
  mgr <- matchitToMatched(m, ranges = gr, keep_mcols = TRUE)
  expect_s4_class(mgr, "MatchedGRanges")
  expect_identical(granges(focal(mgr)), granges(gr)[1:10])
  expect_identical(granges(pool(mgr)), granges(gr)[11:50])
  idx <- indices(mgr)
  expect_identical(sort(idx), sort(pos))
  expect_identical(granges(mgr), granges(gr)[11:50][idx])
  expect_identical(names(mcols(mgr)), c("score", "treat", "x", "g"))
  expect_identical(mgr$score, gr$score[11:50][idx])

  mgr <- matchitToMatched(m, ranges = gr, keep_mcols = FALSE)
  expect_identical(indices(mgr), idx)
  expect_identical(names(mcols(mgr)), c("treat", "x", "g"))
  expect_identical(mgr$x, d$x[11:50][idx])

  gi <- InteractionSet::GInteractions(anchor1 = 1:50, anchor2 = 50:1,
                                      regions = granges(gr))
  mgi <- matchitToMatched(m, ranges = gi, keep_mcols = FALSE)
  expect_s4_class(mgi, "MatchedGInteractions")
  expect_identical(indices(mgi), idx)
  expect_identical(anchors(mgi, "first"), anchors(gi, "first")[11:50][idx])
  expect_identical(anchors(mgi, "second"), anchors(gi, "second")[11:50][idx])
})

test_that("matchitToMatched builds GRanges from range columns in the data", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("mariner")
  d <- matchit_data()
  d$seqnames <- "chr3"
  d$start <- 1:50 * 10L
  d$end <- d$start + 4L
  m <- MatchIt::matchit(treat ~ x + g, data = d,
                        method = "nearest", replace = FALSE)
  pos <- as.integer(m$match.matrix[, 1]) - 10L

  mgr <- matchitToMatched(m)
  expect_s4_class(mgr, "MatchedGRanges")
  idx <- indices(mgr)
  expect_identical(sort(idx), sort(pos))
  expect_identical(start(focal(mgr)), d$start[1:10])
  expect_identical(start(pool(mgr)), d$start[11:50])
  expect_identical(start(mgr), d$start[11:50][idx])
  expect_identical(as.character(seqnames(mgr)), rep("chr3", 10))
})
