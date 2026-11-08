library(nullranges)
library(GenomicRanges)

# Each test replays the random draws bootRanges makes (same calls, same
# order) and moves the features by brute force, so the expected ranges
# come from the algorithm, not from a recording of the output.

boot_df <- function(br) {
  data.frame(seqnames = as.character(seqnames(br)),
             start = start(br), end = end(br),
             strand = as.character(strand(br)),
             id = br$id,
             iter = as.integer(as.character(br$iter)))
}

# expected data frames are built from doubles; match boot_df's types
finish <- function(df) {
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  rownames(df) <- NULL
  df
}

# Every feature of y overlapping bait block i is shifted by
# target_start[i] - bait_start[i], put on target_chr[i] and clipped to that
# chromosome; a feature left with no bases is dropped.
move_features <- function(y, bait_chr, bait_start, bait_width,
                          target_chr, target_start) {
  ychr <- as.character(seqnames(y))
  out <- lapply(seq_along(bait_start), function(i) {
    bait_end <- bait_start[i] + bait_width[i] - 1
    hit <- which(ychr == bait_chr[i] &
                 start(y) <= bait_end & end(y) >= bait_start[i])
    d <- target_start[i] - bait_start[i]
    s <- pmax(start(y)[hit] + d, 1)
    e <- pmin(end(y)[hit] + d, seqlengths(y)[[target_chr[i]]])
    keep <- s <= e
    data.frame(seqnames = rep(target_chr[i], sum(keep)),
               start = s[keep], end = e[keep],
               strand = as.character(strand(y))[hit][keep],
               id = y$id[hit][keep])
  })
  do.call(rbind, out)
}

# Tiles of width L_b along each chromosome, the last one cut at the end.
tile_starts <- function(lens, L_b) {
  n <- ceiling(lens / L_b)
  data.frame(chr = rep(names(lens), n),
             start = unlist(lapply(n, function(k) (seq_len(k) - 1) * L_b + 1),
                            use.names = FALSE))
}

make_y <- function(lens) {
  y <- GRanges(rep(c("chr1", "chr2"), c(8, 5)),
               IRanges(c(5, 40, 95, 150, 210, 290, 350, 385,
                         10, 60, 120, 180, 215),
                       width = c(10, 20, 25, 10, 30, 20, 10, 15,
                                 15, 10, 30, 15, 10)),
               strand = rep(c("+", "-", "*"), length.out = 13),
               seqlengths = lens)
  # across-chromosome bootstrapping requires sort(y) order
  y <- sort(y)
  y$id <- seq_along(y)
  y
}

test_that("bootstrap within chromosome matches a replay of its block draws", {
  # chromosome lengths are multiples of blockLength, so no tile runs past
  # a chromosome end
  y <- make_y(c(chr1 = 400, chr2 = 300))
  L_b <- 100
  R <- 3
  set.seed(11)
  br <- bootRanges(y, blockLength = L_b, R = R,
                   type = "bootstrap", withinChrom = TRUE)

  set.seed(11)
  expected <- do.call(rbind, lapply(seq_len(R), function(r) {
    per_chr <- lapply(seqlevels(y), function(chr) {
      n <- seqlengths(y)[[chr]] / L_b
      bait <- round(runif(n, 1, (n - 1) * L_b + 1))
      # bait block i lands on the i-th tile of the same chromosome
      move_features(y, rep(chr, n), bait, rep(L_b, n),
                    rep(chr, n), (seq_len(n) - 1) * L_b + 1)
    })
    cbind(do.call(rbind, per_chr), iter = r)
  }))
  expect_identical(boot_df(br), finish(expected))
})

test_that("permute within chromosome matches a replay of sample(n)", {
  y <- make_y(c(chr1 = 400, chr2 = 300))
  L_b <- 100
  R <- 3
  set.seed(12)
  br <- bootRanges(y, blockLength = L_b, R = R,
                   type = "permute", withinChrom = TRUE)

  set.seed(12)
  expected <- do.call(rbind, lapply(seq_len(R), function(r) {
    per_chr <- lapply(seqlevels(y), function(chr) {
      n <- seqlengths(y)[[chr]] / L_b
      perm <- sample(n)
      r_chr <- y[seqnames(y) == chr]
      # a feature moves with the tile holding its start: tile b -> perm[b]
      b <- (start(r_chr) - 1) %/% L_b + 1
      d <- (perm[b] - b) * L_b
      data.frame(seqnames = chr,
                 start = pmax(start(r_chr) + d, 1),
                 end = pmin(end(r_chr) + d, seqlengths(y)[[chr]]),
                 strand = as.character(strand(r_chr)),
                 id = r_chr$id)
    })
    cbind(do.call(rbind, per_chr), iter = r)
  }))
  expect_identical(boot_df(br), finish(expected))
})

test_that("bootstrap across chromosomes matches a replay of its block draws", {
  # cut last tiles (401-450, 201-230) make some moved features clip or drop
  y <- make_y(c(chr1 = 450, chr2 = 230))
  L_b <- 100
  R <- 3
  lens <- seqlengths(y)
  tiles <- tile_starts(lens, L_b)
  set.seed(13)
  br <- bootRanges(y, blockLength = L_b, R = R,
                   type = "bootstrap", withinChrom = FALSE)

  set.seed(13)
  expected <- do.call(rbind, lapply(seq_len(R), function(r) {
    n_per_chr <- ceiling(lens / L_b)
    n <- sum(n_per_chr)
    bait_chr <- sample(names(lens), n, replace = TRUE, prob = lens)
    bait <- round(runif(n, 1, (n_per_chr[bait_chr] - 1) * L_b + 1))
    # bait block i lands on tile i of the genome-wide tiling
    cbind(move_features(y, bait_chr, bait, rep(L_b, n),
                        tiles$chr, tiles$start),
          iter = r)
  }))
  expect_identical(boot_df(br), finish(expected))
  # the seed exercises clipping at chromosome ends
  expect_true(any(expected$end - expected$start + 1 != width(y)[expected$id]))
})

test_that("permute across chromosomes matches a replay of sample(n)", {
  y <- make_y(c(chr1 = 450, chr2 = 230))
  L_b <- 100
  R <- 3
  lens <- seqlengths(y)
  tiles <- tile_starts(lens, L_b)
  set.seed(14)
  br <- bootRanges(y, blockLength = L_b, R = R,
                   type = "permute", withinChrom = FALSE)

  ychr <- as.character(seqnames(y))
  b <- vapply(seq_along(y), function(j) {
    which(tiles$chr == ychr[j] & tiles$start <= start(y)[j] &
          tiles$start + L_b > start(y)[j])
  }, integer(1))
  set.seed(14)
  expected <- do.call(rbind, lapply(seq_len(R), function(r) {
    perm <- sample(nrow(tiles))
    to <- perm[b]
    d <- tiles$start[to] - tiles$start[b]
    s <- pmax(start(y) + d, 1)
    e <- pmin(end(y) + d, lens[tiles$chr[to]])
    keep <- s <= e
    data.frame(seqnames = tiles$chr[to][keep], start = s[keep],
               end = unname(e[keep]),
               strand = as.character(strand(y))[keep],
               id = y$id[keep], iter = r)
  }))
  expect_identical(boot_df(br), finish(expected))
  # each feature is used at most once per iteration
  expect_false(anyDuplicated(paste(br$iter, br$id)) > 0)
})

seg_y <- function() {
  y <- GRanges(rep(c("chr1", "chr2"), c(9, 6)),
               IRanges(c(10, 55, 130, 190, 260, 330, 390, 470, 560,
                         20, 90, 150, 240, 300, 370),
                       width = c(12, 20, 15, 30, 10, 25, 20, 15, 30,
                                 10, 20, 15, 25, 10, 20)),
               strand = rep(c("+", "-"), length.out = 15),
               seqlengths = c(chr1 = 600, chr2 = 400))
  y$id <- seq_along(y)
  y
}
seg_2state <- function() {
  GRanges(c("chr1", "chr1", "chr2", "chr2"),
          IRanges(c(1, 404, 1, 201), c(403, 600, 200, 400)),
          state = c(1, 2, 2, 1))
}

test_that("segmented bootstrap with proportional blocks matches a replay", {
  y <- seg_y()
  seg <- seg_2state()
  L_b <- 100
  R <- 2
  set.seed(15)
  br <- bootRanges(y, blockLength = L_b, R = R, seg = seg)

  # state 1 covers 603 of 1000 bp -> blocks of round(60.3) = 60;
  # state 2 covers 397 -> round(39.7) = 40
  L_bs <- c(60, 40)
  seg_chr <- as.character(seqnames(seg))
  set.seed(15)
  expected <- do.call(rbind, lapply(seq_len(R), function(r) {
    per_state <- lapply(1:2, function(s) {
      i <- which(seg$state == s)
      w <- width(seg)[i]
      nb <- ceiling(w / L_bs[s])
      pick <- sample(length(i), sum(nb), replace = TRUE, prob = w)
      rmin <- start(seg)[i][pick]
      bait <- round(runif(sum(nb), rmin, rmin + (nb[pick] - 1) * L_bs[s]))
      # targets tile each segment of this state from its start
      target <- unlist(lapply(seq_along(i), function(k) {
        start(seg)[i][k] + (seq_len(nb[k]) - 1) * L_bs[s]
      }))
      list(bait_chr = seg_chr[i][pick], bait = bait,
           width = rep(L_bs[s], sum(nb)),
           target_chr = rep(seg_chr[i], nb), target = target)
    })
    blocks <- lapply(names(per_state[[1]]), function(nm) {
      unlist(lapply(per_state, `[[`, nm))
    })
    names(blocks) <- names(per_state[[1]])
    cbind(move_features(y, blocks$bait_chr, blocks$bait, blocks$width,
                        blocks$target_chr, blocks$target),
          iter = r)
  }))
  expect_identical(boot_df(br), finish(expected))
})

test_that("segmented bootstrap without proportional blocks matches a replay", {
  y <- seg_y()
  seg <- seg_2state()
  L_b <- 100
  R <- 4
  set.seed(16)
  br <- bootRanges(y, blockLength = L_b, R = R, seg = seg,
                   proportionLength = FALSE)

  seg_chr <- as.character(seqnames(seg))
  nb <- ceiling(width(seg) / L_b)
  targets <- lapply(seq_along(seg), function(k) {
    start(seg)[k] + (seq_len(nb[k]) - 1) * L_b
  })
  branch <- character(0)
  set.seed(16)
  expected <- do.call(rbind, lapply(seq_len(R), function(r) {
    pick <- sample(length(seg), sum(nb), replace = TRUE, prob = width(seg))
    rmin <- start(seg)[pick]
    bait <- round(runif(sum(nb), rmin, rmin + (nb[pick] - 1) * L_b))
    per_state <- lapply(1:2, function(s) {
      drawn <- which(seg$state[pick] == s)
      n_target <- sum(nb[seg$state == s])
      # too few draws for this state: resample them; too many: keep the first
      if (length(drawn) < n_target) {
        branch <<- c(branch, "up")
        drawn <- drawn[sample(length(drawn), n_target, replace = TRUE)]
      } else {
        branch <<- c(branch, "down")
        drawn <- drawn[seq_len(n_target)]
      }
      list(bait_chr = seg_chr[pick][drawn], bait = bait[drawn],
           target_chr = rep(seg_chr[seg$state == s], nb[seg$state == s]),
           target = unlist(targets[seg$state == s]))
    })
    blocks <- lapply(names(per_state[[1]]), function(nm) {
      unlist(lapply(per_state, `[[`, nm))
    })
    names(blocks) <- names(per_state[[1]])
    cbind(move_features(y, blocks$bait_chr, blocks$bait,
                        rep(L_b, length(blocks$bait)),
                        blocks$target_chr, blocks$target),
          iter = r)
  }))
  expect_identical(boot_df(br), finish(expected))
  # the seed exercises both the up- and the down-sampling branch
  expect_setequal(branch, c("up", "down"))
})
