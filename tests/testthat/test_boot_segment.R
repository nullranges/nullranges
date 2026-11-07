library(nullranges)
library(GenomicRanges)

test_that("oneRegionSegment splits a chromosome around one region", {
  x <- GRanges("chr2", IRanges(101, 150))
  seg <- oneRegionSegment(x, seqlength = 1000)
  expect_identical(as.character(seqnames(seg)), rep("chr2", 3))
  expect_identical(start(seg), c(1L, 101L, 151L))
  expect_identical(end(seg), c(100L, 150L, 1000L))
  expect_equal(seg$state, c(1, 2, 1))

  expect_error(oneRegionSegment(x), "seqlength")
  expect_error(oneRegionSegment(c(x, x), seqlength = 1000))
})

test_that("reduceSegment merges touching ranges of the same state", {
  x <- GRanges(c("chr2", "chr1", "chr1", "chr1", "chr1", "chr2"),
               IRanges(c(1, 1, 11, 21, 31, 51), c(50, 10, 20, 30, 40, 60)),
               score = 1:6)
  seqlevels(x) <- c("chr2", "chr1")
  # by hand: state 1 is chr1 1-10 + 11-20 (touching -> 1-20), chr1 31-40
  # and chr2 1-50; state 2 is chr1 21-30 and chr2 51-60
  expected_chr <- c("chr1", "chr1", "chr1", "chr2", "chr2")
  expected_start <- c(1L, 21L, 31L, 1L, 51L)
  expected_end <- c(20L, 30L, 40L, 50L, 60L)
  expected_state <- c(1L, 2L, 1L, 1L, 2L)

  x$state <- c(1, 1, 1, 2, 1, 2)
  seg <- reduceSegment(x)
  expect_identical(seqlevels(seg), c("chr1", "chr2"))
  expect_identical(as.character(seqnames(seg)), expected_chr)
  expect_identical(start(seg), expected_start)
  expect_identical(end(seg), expected_end)
  expect_identical(seg$state, expected_state)
  expect_identical(names(mcols(seg)), "state")

  # the same states under another column name
  x$state <- NULL
  x$label <- c(1, 1, 1, 2, 1, 2)
  seg2 <- reduceSegment(x, col = "label")
  expect_identical(seg2, seg)
})

# Two density levels, 30 and 5 features per 200 bp, in 4 kb regions:
# chr1 dense / sparse / dense, chr2 sparse / dense, nothing on chr3.
density_x <- function() {
  set.seed(30)
  region <- data.frame(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"),
                       start = c(1, 4001, 8001, 1, 4001),
                       per_tile = c(30, 5, 30, 5, 30))
  # features of width 10 lie wholly inside their region
  starts <- unlist(lapply(seq_len(nrow(region)), function(i) {
    round(runif(20 * region$per_tile[i], region$start[i],
                region$start[i] + 3990))
  }))
  x <- GRanges(rep(region$chr, 20 * region$per_tile),
               IRanges(starts, width = 10),
               seqlengths = c(chr1 = 12000, chr2 = 8000, chr3 = 5000))
  sort(x)
}

# counts column by brute force: features overlapping a segment by >= 8 bp,
# per L_s bp of segment
brute_counts <- function(seg, x, L_s) {
  xchr <- as.character(seqnames(x))
  vapply(seq_along(seg), function(i) {
    same <- xchr == as.character(seqnames(seg))[i]
    ov <- pmin(end(x)[same], end(seg)[i]) -
      pmax(start(x)[same], start(seg)[i]) + 1
    sum(ov >= 8) / width(seg)[i] * L_s
  }, numeric(1))
}

expect_density_segments <- function(seg, x, chr, start, end, dense, L_s) {
  expect_identical(as.character(seqnames(seg)), chr)
  expect_identical(start(seg), as.integer(start))
  expect_identical(end(seg), as.integer(end))
  # state labels are arbitrary: dense and sparse ranges get one each
  expect_true(all(seg$state %in% 1:2))
  expect_length(unique(seg$state[dense]), 1)
  expect_length(unique(seg$state[!dense]), 1)
  expect_false(seg$state[dense][1] == seg$state[!dense][1])
  expect_equal(seg$counts, brute_counts(seg, x, L_s))
}

test_that("segmentDensity with CBS separates dense from sparse regions", {
  skip_if_not_installed("DNAcopy")
  x <- density_x()
  set.seed(31)
  capture.output(seg <- segmentDensity(x, n = 2, L_s = 200, type = "cbs"))
  expect_density_segments(seg, x,
                          chr = rep(c("chr1", "chr2"), c(3, 2)),
                          start = c(1, 4001, 8001, 1, 4001),
                          end = c(4000, 8000, 12000, 4000, 8000),
                          dense = c(TRUE, FALSE, TRUE, FALSE, TRUE),
                          L_s = 200)

  # tiles and counts do not depend on input order, so only a warning differs
  on_chr1 <- which(as.character(seqnames(x)) == "chr1")
  x_unsorted <- x[c(rev(on_chr1), seq_along(x)[-on_chr1])]
  set.seed(31)
  expect_warning(
    capture.output(seg2 <- segmentDensity(x_unsorted, n = 2, L_s = 200,
                                          type = "cbs")),
    "unsorted"
  )
  expect_identical(seg2, seg)
})

test_that("segmentDensity with CBS leaves out excluded bases", {
  skip_if_not_installed("DNAcopy")
  x <- density_x()
  exclude <- GRanges(c("chr1", "chr2"), IRanges(c(5001, 6101), c(5199, 6300)),
                     seqlengths = seqlengths(x))
  set.seed(32)
  capture.output(seg <- segmentDensity(x, n = 2, L_s = 200,
                                       exclude = exclude, type = "cbs"))
  # chr1 tile 5001-5200 keeps one base (5200), under L_s / 100 = 2 bp, so it
  # is dropped; chr2 tiles 6001-6200 and 6201-6400 keep 6001-6100, 6301-6400
  expect_density_segments(seg, x,
                          chr = rep(c("chr1", "chr2"), c(4, 3)),
                          start = c(1, 4001, 5201, 8001, 1, 4001, 6301),
                          end = c(4000, 5000, 8000, 12000, 4000, 6100, 8000),
                          dense = c(TRUE, FALSE, FALSE, TRUE,
                                    FALSE, TRUE, TRUE),
                          L_s = 200)
  expect_false(any(overlapsAny(seg, exclude)))
})

test_that("segmentDensity with an HMM separates dense from sparse regions", {
  skip_if_not_installed("RcppHMM")
  x <- density_x()
  # initGHMM draws EM's start from R's RNG. 7 of seeds 1-100 stall with two
  # equal states; seed 1's start converges in 5 iterations, and 800 starts
  # perturbed around it (noise sd up to 0.3) give the same segments, so
  # platform rounding in RcppHMM cannot change them.
  set.seed(1)
  capture.output(seg <- segmentDensity(x, n = 2, L_s = 200, type = "hmm"))
  expect_density_segments(seg, x,
                          chr = rep(c("chr1", "chr2"), c(3, 2)),
                          start = c(1, 4001, 8001, 1, 4001),
                          end = c(4000, 8000, 12000, 4000, 8000),
                          dense = c(TRUE, FALSE, TRUE, FALSE, TRUE),
                          L_s = 200)
})
