library(nullranges)
library(GenomicRanges)

prop_y <- function() {
  set.seed(1)
  lens <- c(chr1 = 100000, chr2 = 61234, chr3 = 35500)
  n <- c(1000, 600, 400)
  starts <- unlist(lapply(1:3, function(i) {
    round(runif(n[i], 1, lens[i] - 200))
  }))
  y <- GRanges(rep(names(lens), n),
               IRanges(starts, width = round(runif(sum(n), 20, 200))),
               strand = sample(c("+", "-", "*"), sum(n), replace = TRUE),
               seqlengths = lens)
  y <- sort(y)
  y$id <- seq_along(y)
  y$score <- runif(length(y))
  y
}
prop_seg <- function() {
  GRanges(c("chr1", "chr1", "chr1", "chr2", "chr2", "chr3"),
          IRanges(c(1, 40001, 70001, 1, 30001, 1),
                  c(40000, 70000, 100000, 30000, 61234, 35500)),
          state = c(1, 2, 1, 2, 1, 2))
}

test_that("every bootstrap variant returns well-formed ranges", {
  y <- prop_y()
  seg <- prop_seg()
  L_b <- 5000
  R <- 5
  lens <- seqlengths(y)
  variants <- list(
    across_boot = list(type = "bootstrap", withinChrom = FALSE),
    across_perm = list(type = "permute", withinChrom = FALSE),
    within_boot = list(type = "bootstrap", withinChrom = TRUE),
    within_perm = list(type = "permute", withinChrom = TRUE),
    seg_prop = list(seg = seg),
    seg_noprop = list(seg = seg, proportionLength = FALSE)
  )
  for (v in names(variants)) {
    args <- c(list(y = y, blockLength = L_b, R = R), variants[[v]])
    set.seed(2)
    br <- do.call(bootRanges, args)
    set.seed(2)
    br2 <- do.call(bootRanges, args)
    info <- paste("variant", v)

    expect_identical(br, br2, info = info)
    expect_s4_class(br, "BootRanges")
    expect_true(is(br, "GRanges"), info = info)
    expect_identical(seqinfo(br), seqinfo(y), info = info)
    expect_identical(names(mcols(br)), c("id", "score", "iter"), info = info)
    expect_s4_class(br$iter, "Rle")
    expect_true(is.factor(runValue(br$iter)), info = info)
    expect_identical(levels(br$iter), as.character(seq_len(R)), info = info)
    expect_false(is.unsorted(as.integer(br$iter)), info = info)
    expect_true(all(tabulate(as.integer(br$iter), R) > 0), info = info)
    # metadata travels with each range
    expect_identical(br$score, y$score[br$id], info = info)
    # within the chromosome
    chr_len <- lens[as.character(seqnames(br))]
    expect_true(all(start(br) >= 1 & end(br) <= chr_len), info = info)
    # block shifts keep widths, except where trimmed at a chromosome end
    kept <- width(br) == width(y)[br$id]
    expect_true(all(kept | start(br) == 1 | end(br) == chr_len), info = info)
    expect_gt(mean(kept), 0.9)
    # and the ranges do move
    expect_gt(mean(start(br) != start(y)[br$id]), 0.5)
  }
})

test_that("within-chromosome variants keep features on their chromosome", {
  y <- prop_y()
  R <- 4
  for (type in c("bootstrap", "permute")) {
    set.seed(3)
    br <- bootRanges(y, blockLength = 5000, R = R, type = type,
                     withinChrom = TRUE)
    expect_identical(as.character(seqnames(br)),
                     as.character(seqnames(y))[br$id], info = type)
  }
  # permuting (br from the last pass) moves each feature at most once per
  # iteration, and exactly once on chr1, whose length is a multiple of
  # blockLength (on chr2 and chr3 a feature moved past the end may be lost)
  on_chr1 <- y$id[as.character(seqnames(y)) == "chr1"]
  for (r in seq_len(R)) {
    ids <- br$id[as.integer(br$iter) == r]
    expect_false(anyDuplicated(ids) > 0)
    expect_true(all(on_chr1 %in% ids))
  }
})

test_that("within-chromosome variants do not need sorted input", {
  y <- prop_y()
  y_shuffled <- y[c(seq(2, length(y), by = 2), seq(1, length(y), by = 2))]
  # blocks are drawn per chromosome in seqlevels order, so shuffling the
  # input only reorders the output within each iteration
  canonical <- function(br) {
    df <- data.frame(iter = as.integer(br$iter), id = br$id,
                     seqnames = as.character(seqnames(br)),
                     start = start(br), end = end(br),
                     strand = as.character(strand(br)))
    df <- df[do.call(order, df), ]
    rownames(df) <- NULL
    df
  }
  for (type in c("bootstrap", "permute")) {
    set.seed(4)
    br <- bootRanges(y_shuffled, blockLength = 5000, R = 2,
                     type = type, withinChrom = TRUE)
    set.seed(4)
    br_sorted <- bootRanges(y, blockLength = 5000, R = 2,
                            type = type, withinChrom = TRUE)
    expect_identical(canonical(br), canonical(br_sorted), info = type)
    expect_false(identical(br$id, br_sorted$id), info = type)
  }
  expect_error(bootRanges(y_shuffled, blockLength = 5000))
})

test_that("permuting across chromosomes uses each feature at most once", {
  y <- prop_y()
  set.seed(5)
  br <- bootRanges(y, blockLength = 5000, R = 4, type = "permute")
  expect_false(anyDuplicated(paste(br$iter, br$id)) > 0)
  # features only go missing when moved into a cut last tile and trimmed away
  per_iter <- tabulate(as.integer(br$iter), 4)
  expect_true(all(per_iter <= length(y)))
  expect_gt(min(per_iter), 0.9 * length(y))
})

test_that("segmented bootstrap ignores type and withinChrom", {
  y <- prop_y()
  seg <- prop_seg()
  set.seed(6)
  br1 <- bootRanges(y, blockLength = 5000, R = 2, seg = seg)
  set.seed(6)
  br2 <- bootRanges(y, blockLength = 5000, R = 2, seg = seg,
                    type = "permute", withinChrom = TRUE)
  expect_identical(br1, br2)
})

test_that("segmented bootstrap mostly keeps features in their state", {
  y <- prop_y()
  seg <- prop_seg()
  state_of <- function(gr) {
    seg$state[findOverlaps(gr, seg, select = "first")]
  }
  for (prop in c(TRUE, FALSE)) {
    set.seed(7)
    br <- bootRanges(y, blockLength = 5000, R = 2, seg = seg,
                     proportionLength = prop)
    same <- state_of(br) == state_of(y)[br$id]
    expect_gt(mean(same), 0.9)
  }
})

test_that("storeBlockLength records the block length as an integer Rle", {
  y <- prop_y()
  set.seed(8)
  br <- bootRanges(y, blockLength = 5000, R = 2, storeBlockLength = TRUE)
  expect_s4_class(br$blockLength, "Rle")
  expect_identical(as.vector(br$blockLength), rep(5000L, length(br)))
  # storing the block length does not change the ranges
  set.seed(8)
  br0 <- bootRanges(y, blockLength = 5000, R = 2)
  expect_null(br0$blockLength)
  expect_identical(granges(br), granges(br0))
})

test_that("iterations with no ranges keep their iter level", {
  # one feature in tile 1 of 2; permuting moves it to tile 2 (excluded)
  # whenever sample(2) starts with 2
  y <- GRanges("chr1", IRanges(10, 19), seqlengths = c(chr1 = 200))
  exclude <- GRanges("chr1", IRanges(101, 200))
  R <- 8
  set.seed(10)
  kept <- vapply(seq_len(R), function(r) sample(2)[1] == 1, logical(1))
  expect_true(any(kept) && !all(kept))
  set.seed(10)
  br <- bootRanges(y, blockLength = 100, R = R, exclude = exclude,
                   type = "permute", withinChrom = TRUE)
  expect_identical(levels(br$iter), as.character(seq_len(R)))
  expect_identical(tabulate(as.integer(br$iter), R), as.integer(kept))
  expect_identical(start(br), rep(10L, sum(kept)))
})

test_that("iterations of exactly one range are labelled 1..R", {
  # a single tile: permuting cannot move the feature
  y <- GRanges("chr1", IRanges(5, 10), seqlengths = c(chr1 = 100))
  br <- bootRanges(y, blockLength = 100, R = 3, type = "permute")
  expect_length(br, 3)
  expect_identical(start(br), rep(5L, 3))
  expect_identical(as.integer(br$iter), 1:3)
  br1 <- bootRanges(y, blockLength = 100)
  expect_identical(levels(br1$iter), "1")
  expect_length(br1, 1)
})

test_that("bootRanges checks its inputs", {
  y <- GRanges("chr1", IRanges(c(5, 50), width = 10),
               seqlengths = c(chr1 = 100))
  y_nolen <- y
  seqlengths(y_nolen) <- NA
  # each call is valid apart from the one input named in its comment
  # no seqlengths
  expect_error(bootRanges(y_nolen, blockLength = 20))
  # block longer than the chromosome
  expect_error(bootRanges(y, blockLength = 101))
  # unsorted input when bootstrapping across chromosomes
  expect_error(bootRanges(rev(y), blockLength = 20))
  # unknown type
  expect_error(bootRanges(y, blockLength = 20, type = "shuffle"))
  # segmentation without a state column
  seg_nostate <- GRanges("chr1", IRanges(1, 100))
  expect_error(bootRanges(y, blockLength = 20, seg = seg_nostate))
})
