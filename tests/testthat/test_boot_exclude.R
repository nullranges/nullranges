library(nullranges)
library(GenomicRanges)

# Excluding happens after the draws of each iteration and draws nothing
# itself, so under the same seed the output with `exclude` must equal the
# output without it, filtered (drop) or cut (trim) by brute force.

excl_df <- function(gr) {
  data.frame(seqnames = as.character(seqnames(gr)),
             start = start(gr), end = end(gr),
             strand = as.character(strand(gr)),
             id = gr$id,
             iter = as.integer(gr$iter))
}

excl_inputs <- function() {
  set.seed(20)
  lens <- c(chr1 = 2000, chr2 = 1000)
  n <- c(40, 20)
  starts <- c(round(runif(n[1], 1, 1940)), round(runif(n[2], 1, 940)))
  y <- GRanges(rep(names(lens), n),
               IRanges(starts, width = round(runif(sum(n), 10, 60))),
               strand = rep(c("+", "-"), length.out = sum(n)),
               seqlengths = lens)
  y <- sort(y)
  y$id <- seq_along(y)
  # short excludes so that some ranges span one and are split in two
  exclude <- GRanges(c("chr1", "chr1", "chr1", "chr2", "chr2"),
                     IRanges(c(301, 1201, 1601, 151, 501),
                             c(380, 1205, 1603, 153, 700)),
                     seqlengths = lens)
  seg <- GRanges(c("chr1", "chr1", "chr2"),
                 IRanges(c(1, 1001, 1), c(1000, 2000, 1000)),
                 state = c(1, 2, 1))
  list(y = y, exclude = exclude, seg = seg)
}

# the three code paths that produce y_prime before excluding
excl_variants <- function(inp) {
  list(
    within_boot = list(type = "bootstrap", withinChrom = TRUE),
    across_perm = list(type = "permute", withinChrom = FALSE),
    seg_prop = list(seg = inp$seg)
  )
}

boot_with <- function(inp, args, seed, ...) {
  set.seed(seed)
  do.call(bootRanges, c(list(y = inp$y, blockLength = 100, R = 3),
                        args, list(...)))
}

test_that("excludeOption = 'drop' removes exactly the ranges hitting exclude", {
  inp <- excl_inputs()
  ex <- inp$exclude
  ex_chr <- as.character(seqnames(ex))
  for (v in names(excl_variants(inp))) {
    args <- excl_variants(inp)[[v]]
    full <- excl_df(boot_with(inp, args, 21))
    br <- boot_with(inp, args, 21, exclude = ex, excludeOption = "drop")

    hits <- vapply(seq_len(nrow(full)), function(i) {
      any(ex_chr == full$seqnames[i] &
          start(ex) <= full$end[i] & end(ex) >= full$start[i])
    }, logical(1))
    expected <- full[!hits, ]
    rownames(expected) <- NULL
    expect_identical(excl_df(br), expected, info = v)
    expect_true(any(hits), info = v)
    expect_false(any(overlapsAny(br, ex)), info = v)
    expect_identical(levels(br$iter), c("1", "2", "3"), info = v)
  }
  # drop is the default
  args <- excl_variants(inp)$across_perm
  expect_identical(boot_with(inp, args, 21, exclude = ex),
                   boot_with(inp, args, 21, exclude = ex,
                             excludeOption = "drop"))
})

test_that("excludeOption = 'trim' removes exactly the excluded bases", {
  inp <- excl_inputs()
  ex <- inp$exclude
  ex_chr <- as.character(seqnames(ex))
  for (v in names(excl_variants(inp))) {
    args <- excl_variants(inp)[[v]]
    full <- excl_df(boot_with(inp, args, 22))
    br <- boot_with(inp, args, 22, exclude = ex, excludeOption = "trim")

    # each range minus its excluded bases, one row per remaining run of bases
    pieces <- lapply(seq_len(nrow(full)), function(i) {
      pos <- full$start[i]:full$end[i]
      on_chr <- which(ex_chr == full$seqnames[i])
      for (k in on_chr) {
        pos <- pos[pos < start(ex)[k] | pos > end(ex)[k]]
      }
      if (length(pos) == 0) return(NULL)
      run <- cumsum(c(1, diff(pos) != 1))
      data.frame(seqnames = full$seqnames[i],
                 start = as.integer(tapply(pos, run, min)),
                 end = as.integer(tapply(pos, run, max)),
                 strand = full$strand[i], id = full$id[i],
                 iter = full$iter[i])
    })
    expected <- do.call(rbind, pieces)
    rownames(expected) <- NULL
    expect_identical(excl_df(br), expected, info = v)
    expect_false(any(overlapsAny(br, ex)), info = v)
    # the seed exercises removal, splitting around an exclude and clipping
    n_pieces <- vapply(pieces, function(p) if (is.null(p)) 0L else nrow(p),
                       integer(1))
    remaining <- vapply(pieces, function(p) {
      if (is.null(p)) 0L else sum(p$end - p$start + 1L)
    }, integer(1))
    clipped <- n_pieces == 1 & remaining < full$end - full$start + 1
    expect_true(any(n_pieces == 0) && any(n_pieces == 2) && any(clipped),
                info = v)
  }
})

test_that("trim needs an unstranded exclude and excludeOption is checked", {
  inp <- excl_inputs()
  ex_plus <- inp$exclude
  strand(ex_plus) <- "+"
  expect_error(bootRanges(inp$y, blockLength = 100, exclude = ex_plus,
                          excludeOption = "trim"))
  expect_error(bootRanges(inp$y, blockLength = 100, exclude = inp$exclude,
                          excludeOption = "cut"))
})
