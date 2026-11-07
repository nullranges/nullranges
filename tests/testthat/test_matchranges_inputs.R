library(nullranges)

## 12 focal and 60 pool rows with numeric, character and logical covariates
inputs_data <- function() {
  set.seed(1)
  n <- 72
  data.frame(is_focal = rep(c(TRUE, FALSE), c(12, 60)),
             x = c(rnorm(12, 1), rnorm(60)),
             g = sample(c("u", "v", "w"), n, replace = TRUE),
             flag = rep(c(TRUE, FALSE), length.out = n))
}

test_that("matchRanges returns the Matched subclass for each input class", {
  df <- inputs_data()
  gr <- GRanges("chr1", IRanges(start = seq_len(nrow(df)) * 10, width = 5))
  mcols(gr) <- df
  gi <- InteractionSet::GInteractions(anchor1 = seq_len(nrow(df)),
                                      anchor2 = rev(seq_len(nrow(df))),
                                      regions = gr)
  mcols(gi) <- df
  inputs <- list(data.frame = df,
                 data.table = data.table::as.data.table(df),
                 DataFrame = DataFrame(df),
                 GRanges = gr,
                 GInteractions = gi)
  expected_class <- c(data.frame = "MatchedDataFrame",
                      data.table = "MatchedDataFrame",
                      DataFrame = "MatchedDataFrame",
                      GRanges = "MatchedGRanges",
                      GInteractions = "MatchedGInteractions")

  for (meth in c("nearest", "stratified")) {
    set.seed(2)
    ref <- suppressMessages(matchRanges(df[df$is_focal, ], df[!df$is_focal, ],
                                        covar = ~x + g, method = meth,
                                        replace = TRUE))
    for (cl in names(inputs)) {
      obj <- inputs[[cl]]
      focal <- obj[obj$is_focal, ]
      pool <- obj[!obj$is_focal, ]
      set.seed(2)
      m <- suppressMessages(matchRanges(focal, pool, covar = ~x + g,
                                        method = meth, replace = TRUE))
      expect_s4_class(m, expected_class[[cl]])
      expect_identical(method(m), meth)
      expect_true(withReplacement(m))
      expect_identical(covariates(m), c("x", "g"))

      ## the container does not change the match
      idx <- indices(m)
      expect_identical(idx, indices(ref))
      expect_equal(as.data.frame(matchedData(m)),
                   as.data.frame(matchedData(ref)))

      if (cl %in% c("GRanges", "GInteractions")) {
        expect_identical(focal(m), focal)
        expect_identical(pool(m), pool)
        expect_identical(matched(m), pool[idx])
      } else {
        expect_identical(focal(m), DataFrame(focal))
        expect_identical(pool(m), DataFrame(pool))
      }
      ## the object itself holds the matched pool rows
      expect_identical(NROW(m), 12L)
      expect_identical(m$x, pool$x[idx])
      expect_identical(m$g, pool$g[idx])
    }
  }
})

test_that("matched GRanges and GInteractions keep the pool's coordinates", {
  df <- inputs_data()
  gr <- GRanges("chr1", IRanges(start = seq_len(nrow(df)) * 10, width = 5))
  mcols(gr) <- df
  set.seed(3)
  m <- suppressMessages(matchRanges(gr[gr$is_focal], gr[!gr$is_focal],
                                    covar = ~x, method = "stratified",
                                    replace = FALSE))
  pool <- gr[!gr$is_focal]
  expect_identical(granges(m), granges(pool)[indices(m)])
  expect_identical(start(m), 10L * (12L + indices(m)))
  expect_false(anyDuplicated(indices(m)) > 0)

  gi <- InteractionSet::GInteractions(anchor1 = seq_len(nrow(df)),
                                      anchor2 = rev(seq_len(nrow(df))),
                                      regions = granges(gr))
  mcols(gi) <- df
  set.seed(3)
  m <- matchRanges(gi[gi$is_focal], gi[!gi$is_focal], covar = ~x,
                   method = "nearest", replace = TRUE)
  pool <- gi[!gi$is_focal]
  expect_identical(anchors(m, "first"), anchors(pool, "first")[indices(m)])
  expect_identical(anchors(m, "second"), anchors(pool, "second")[indices(m)])
})

test_that("matchedData stacks the focal, matched, pool and unmatched rows", {
  df <- inputs_data()
  focal <- df[df$is_focal, ]
  pool <- df[!df$is_focal, ]
  set.seed(4)
  m <- matchRanges(focal, pool, covar = ~x + g + flag,
                   method = "nearest", replace = TRUE)
  idx <- indices(m)
  md <- as.data.frame(matchedData(m))
  expect_identical(names(md), c("id", "x", "g", "flag", "ps", "set"))
  expect_identical(covariates(m), c("x", "g", "flag"))

  unmatched <- setdiff(seq_len(60), idx)
  expect_identical(md$set, rep(c("focal", "matched", "pool", "unmatched"),
                               c(12, 12, 60, length(unmatched))))
  expect_identical(as.character(md$id),
                   rep(c("1", "0"), c(12, 72 + length(unmatched))))

  rows <- c(seq_len(12), 12 + idx, 12 + seq_len(60), 12 + unmatched)
  for (v in c("x", "g", "flag")) {
    expect_identical(md[[v]], df[[v]][rows])
  }

  ## scores are the fitted values of the logistic regression on all three
  fit <- glm(is_focal ~ x + g + flag, family = binomial("logit"), data = df)
  expect_equal(md$ps, unname(fitted(fit))[rows], tolerance = 1e-8)
})

test_that("character and factor covariates match alike; logical ones work too", {
  df <- inputs_data()
  focal <- df[df$is_focal, ]
  pool <- df[!df$is_focal, ]

  set.seed(5)
  m_chr <- matchRanges(focal, pool, covar = ~g, method = "nearest")
  ## the level scores u/v/w = 4/24, 5/22 and 3/26 are distinct, so every
  ## focal element matches a pool element of its own level
  expect_identical(pool$g[indices(m_chr)], focal$g)
  focal$g <- factor(focal$g, levels = c("u", "v", "w"))
  pool$g <- factor(pool$g, levels = c("u", "v", "w"))
  set.seed(5)
  m_fct <- matchRanges(focal, pool, covar = ~g, method = "nearest")
  expect_identical(indices(m_chr), indices(m_fct))
  expect_identical(matchedData(m_chr)$ps, matchedData(m_fct)$ps)

  ## a single logical covariate: focal TRUE/FALSE = 8/4 and pool 20/40 score
  ## 8/28 and 4/44. Each focal element draws among the pool elements of its
  ## own level, the lower score first: the 4 FALSE focal from pool 1:40, then
  ## the 8 TRUE focal from pool 41:60
  focal <- data.frame(flag = rep(c(TRUE, FALSE), c(8, 4)))
  pool <- data.frame(flag = rep(c(FALSE, TRUE), c(40, 20)))
  set.seed(6)
  m <- matchRanges(focal, pool, covar = ~flag, method = "nearest")
  md <- matchedData(m)
  expect_equal(md$ps[md$set == "focal"], rep(c(8 / 28, 4 / 44), c(8, 4)),
               tolerance = 1e-6)
  set.seed(6)
  low <- sample(40, 4, replace = TRUE)
  high <- 40L + sample(20, 8, replace = TRUE)
  expect_identical(indices(m), c(high, low))
})

test_that("seeded runs repeat exactly", {
  df <- inputs_data()
  focal <- df[df$is_focal, ]
  pool <- df[!df$is_focal, ]
  for (meth in c("nearest", "stratified")) {
    set.seed(7)
    a <- suppressMessages(matchRanges(focal, pool, ~x + g, method = meth))
    set.seed(7)
    b <- suppressMessages(matchRanges(focal, pool, ~x + g, method = meth))
    expect_identical(indices(a), indices(b))
    expect_equal(as.data.frame(matchedData(a)), as.data.frame(matchedData(b)))
  }
})
