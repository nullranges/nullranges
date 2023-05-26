library(nullranges)
library(MatchIt)

test_that("Coerce MatchIt to Matched subclasses", {
  
  ## Matched and MatchedDataFrame ----------------------------------------------
  ## Create example data.frame dataset
  set.seed(123)
  df <- makeExampleMatchedDataSet(type="data.frame")
  
  ## MatchIt
  set.seed(123)
  m <- matchit(formula=feature1 ~ feature2 + feature3,
               data=df,
               method='nearest',
               replace=FALSE)
  
  expect_s4_class(as_Matched(m), "Matched")
  
  ## MatchedDataFrame
  mdf <- as_MatchedDataFrame(m)
  
  expect_s4_class(mdf, "MatchedDataFrame")
  
  expect_s4_class(focal(mdf), "DataFrame")
  expect_s4_class(pool(mdf), "DataFrame")
  expect_s4_class(matched(mdf), "DataFrame")
  expect_s4_class(unmatched(mdf), "DataFrame")
  
  ## matched accessor pulls delegate which is the
  ## same as "matched" indices from the pool set.
  identical(matched(mdf), pool(mdf)[indices(mdf),]) |>
    expect_true()
  identical(unmatched(mdf),
            pool(mdf)[indices(mdf, set="unmatched"),]) |>
    expect_true()
  identical(focal(mdf), focal(mdf)[indices(mdf, set="focal"),]) |>
    expect_true()
  identical(pool(mdf), pool(mdf)[indices(mdf, set="pool"),]) |>
    expect_true()
  
  
  ## MatchedGRanges ------------------------------------------------------------
  ## Create example GRanges dataset
  set.seed(123)
  gr <- makeExampleMatchedDataSet(type="GRanges")
  
  ## MatchIt
  set.seed(123)
  m <- matchit(formula=feature1 ~ feature2 + feature3,
               data=as.data.frame(gr),
               method='nearest',
               replace=FALSE)
  
  mgr <- as_MatchedGRanges(m, ranges=NULL, keep_mcols=TRUE)
  
  expect_s4_class(mgr, "MatchedGRanges")
  
  expect_s4_class(focal(mgr), "GRanges")
  expect_s4_class(pool(mgr), "GRanges")
  expect_s4_class(matched(mgr), "GRanges")
  expect_s4_class(unmatched(mgr), "GRanges")
  
  ## matched accessor pulls delegate which is the
  ## same as "matched" indices from the pool set.
  identical(matched(mgr), pool(mgr)[indices(mgr),]) |>
    expect_true()
  identical(unmatched(mgr),
            pool(mgr)[indices(mgr, set="unmatched"),]) |>
    expect_true()
  identical(focal(mgr), focal(mgr)[indices(mgr, set="focal"),]) |>
    expect_true()
  identical(pool(mgr), pool(mgr)[indices(mgr, set="pool"),]) |>
    expect_true()

  
  ## MatchIt without ranges
  set.seed(123)
  m <- matchit(formula=feature1 ~ feature2 + feature3,
               data=as.data.frame(mcols(gr)),
               method='nearest',
               replace=FALSE)
  
  mgr <- as_MatchedGRanges(m, ranges=gr, keep_mcols=TRUE)
  expect_identical(ncol(mcols(mgr)), 6L)
  
  mgr <- as_MatchedGRanges(m, ranges=gr, keep_mcols=FALSE)
  
  expect_identical(ncol(mcols(mgr)), 3L)
  expect_s4_class(mgr, "MatchedGRanges")
  
  expect_s4_class(focal(mgr), "GRanges")
  expect_s4_class(pool(mgr), "GRanges")
  expect_s4_class(matched(mgr), "GRanges")
  expect_s4_class(unmatched(mgr), "GRanges")
  
  ## matched accessor pulls delegate which is the
  ## same as "matched" indices from the pool set.
  identical(matched(mgr), pool(mgr)[indices(mgr),]) |>
    expect_true()
  identical(unmatched(mgr),
            pool(mgr)[indices(mgr, set="unmatched"),]) |>
    expect_true()
  identical(focal(mgr), focal(mgr)[indices(mgr, set="focal"),]) |>
    expect_true()
  identical(pool(mgr), pool(mgr)[indices(mgr, set="pool"),]) |>
    expect_true()
  
  
  ## MatchedGInteractions ------------------------------------------------------
  ## Create example GInteractions dataset
  set.seed(123)
  gi <- makeExampleMatchedDataSet(type="GInteractions")
  
  ## MatchIt
  set.seed(123)
  m <- matchit(formula=feature1 ~ feature2 + feature3,
               data=as.data.frame(gi),
               method='nearest',
               replace=FALSE)
  
  mgi <- as_MatchedGInteractions(m, interactions=NULL, keep_mcols=TRUE)
  
  expect_s4_class(mgi, "MatchedGInteractions")
  
  expect_s4_class(focal(mgi), "GInteractions")
  expect_s4_class(pool(mgi), "GInteractions")
  expect_s4_class(matched(mgi), "GInteractions")
  expect_s4_class(unmatched(mgi), "GInteractions")
  
  ## matched accessor pulls delegate which is the
  ## same as "matched" indices from the pool set.
  identical(matched(mgi), pool(mgi)[indices(mgi),]) |>
    expect_true()
  identical(unmatched(mgi),
            pool(mgi)[indices(mgi, set="unmatched"),]) |>
    expect_true()
  identical(focal(mgi), focal(mgi)[indices(mgi, set="focal"),]) |>
    expect_true()
  identical(pool(mgi), pool(mgi)[indices(mgi, set="pool"),]) |>
    expect_true()
  
  
  ## MatchIt without ranges
  set.seed(123)
  m <- matchit(formula=feature1 ~ feature2 + feature3,
               data=as.data.frame(mcols(gi)),
               method='nearest',
               replace=FALSE)
  
  mgi <- as_MatchedGInteractions(m, interactions=gi, keep_mcols=TRUE)
  expect_identical(ncol(mcols(mgi)), 6L)
  
  mgi <- as_MatchedGInteractions(m, interactions=gi, keep_mcols=FALSE)
  
  expect_identical(ncol(mcols(mgi)), 3L)
  expect_s4_class(mgi, "MatchedGInteractions")
  
  expect_s4_class(focal(mgi), "GInteractions")
  expect_s4_class(pool(mgi), "GInteractions")
  expect_s4_class(matched(mgi), "GInteractions")
  expect_s4_class(unmatched(mgi), "GInteractions")
  
  ## matched accessor pulls delegate which is the
  ## same as "matched" indices from the pool set.
  identical(matched(mgi), pool(mgi)[indices(mgi),]) |>
    expect_true()
  identical(unmatched(mgi),
            pool(mgi)[indices(mgi, set="unmatched"),]) |>
    expect_true()
  identical(focal(mgi), focal(mgi)[indices(mgi, set="focal"),]) |>
    expect_true()
  identical(pool(mgi), pool(mgi)[indices(mgi, set="pool"),]) |>
    expect_true()
  
  
  ## Testing the full function -------------------------------------------------
  
  ## Create example GRanges dataset
  set.seed(123)
  x <- makeExampleMatchedDataSet(type="GRanges")
  
  ## Ranges without mcols
  ranges <- x
  mcols(ranges) <- NULL
  
  set.seed(123)
  ans <- as.data.frame(mcols(x)) |>
    matchit(formula=feature1 ~ feature2 + feature3,
            data=_,
            method='nearest',
            replace=FALSE) |>
    matchitToMatched(ranges=ranges)
  
  expect_identical(dim(mcols(ans)), c(500L, 3L))
  expect_s4_class(ans, "MatchedGRanges")
  
  set.seed(123)
  ans <- as.data.frame(x) |>
    matchit(formula=feature1 ~ feature2 + feature3,
            data=_,
            method='nearest',
            replace=FALSE) |>
    matchitToMatched()
  
  expect_identical(dim(mcols(ans)), c(500L, 3L))
  expect_s4_class(ans, "MatchedGRanges")
  
  
  ## Create example GInteractions dataset
  set.seed(123)
  x <- makeExampleMatchedDataSet(type="GInteractions")
  
  ## Interactions without mcols
  interactions <- x
  mcols(interactions) <- NULL
  
  set.seed(123)
  ans <- as.data.frame(mcols(x)) |>
    matchit(formula=feature1 ~ feature2 + feature3,
            data=_,
            method='nearest',
            replace=FALSE) |>
    matchitToMatched(ranges=interactions)
  
  expect_identical(dim(mcols(ans)), c(500L, 3L))
  expect_s4_class(ans, "MatchedGInteractions")
  
  set.seed(123)
  ans <- as.data.frame(x) |>
    matchit(formula=feature1 ~ feature2 + feature3,
            data=_,
            method='nearest',
            replace=FALSE) |>
    matchitToMatched()
  
  expect_identical(dim(mcols(ans)), c(500L, 3L))
  expect_s4_class(ans, "MatchedGInteractions")
  
  
  ## Create example data.frame dataset
  set.seed(123)
  x <- makeExampleMatchedDataSet(type="data.table")
  
  set.seed(123)
  ans <- x |>
    matchit(formula=feature1 ~ feature2 + feature3,
            data=_,
            method='nearest',
            replace=FALSE) |>
    matchitToMatched()
  
  expect_identical(dim(ans), c(500L, 3L))
  expect_s4_class(ans, "MatchedDataFrame")
  
})
