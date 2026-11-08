library(nullranges)
library(GenomicRanges)

## matchRanges() with replace neither TRUE nor FALSE: no matching method
## runs, and the matched set is gathered from a variable none made, inside
## a data.table subset, so data.table's error names the columns of the
## pool's covariate table. Nothing is drawn before it.
test_that("replace not TRUE or FALSE fails gathering the matched set", {
  set.seed(1)
  n <- 72
  df <- data.frame(is_focal = rep(c(TRUE, FALSE), c(12, 60)),
                   x = c(rnorm(12, 1), rnorm(60)),
                   log10_width = runif(n))
  gr <- GRanges("chr1", IRanges(start = seq_len(n) * 10, width = 5))
  mcols(gr) <- df
  gi <- InteractionSet::GInteractions(anchor1 = seq_len(n),
                                      anchor2 = rev(seq_len(n)),
                                      regions = gr)
  mcols(gi) <- df
  inputs <- list(data.frame = df,
                 data.table = data.table::as.data.table(df),
                 DataFrame = DataFrame(df),
                 GRanges = gr,
                 GInteractions = gi)
  expected <- c(x = "Object 'mdt' not found amongst [id, x, ps]",
                `x + log10_width` =
                  "Object 'mdt' not found. Perhaps you intended [log10_width]")
  for (cl in names(inputs)) for (meth in c("nearest", "rejection", "stratified"))
    for (replace in list(NA, logical(0), c(TRUE, FALSE)))
      for (cv in names(expected)) {
        obj <- inputs[[cl]]
        focal <- obj[obj$is_focal, ]
        pool <- obj[!obj$is_focal, ]
        set.seed(2)
        seed <- get(".Random.seed", envir = globalenv())
        err <- tryCatch(
          suppressMessages(matchRanges(focal, pool,
                                       covar = as.formula(paste("~", cv)),
                                       method = meth, replace = replace)),
          error = function(e) e)
        info <- paste(cl, meth, deparse(replace), cv)
        expect_s3_class(err, "error")
        expect_identical(conditionMessage(err), expected[[cv]], info = info)
        ## nothing drawn (compared as text, which testthat reports cleanly)
        expect_identical(sprintf("%d", get(".Random.seed", envir = globalenv())),
                         sprintf("%d", seed), info = info)
      }
})
