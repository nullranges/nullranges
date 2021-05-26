library(nullranges)
test_that("trim denyOption works as expected", {

  x <- GRanges("chr1", IRanges(1 + 2 * 0:45, width=10),
               strand=rep(c("+","-"),length=46),
               seqlengths=c(chr1=100))
  deny <- GRanges("chr1", IRanges(36,45),
                  seqlengths=c(chr1=100))
  x <- sort(x)
  x$block <- 1
  b <- bootRanges(x, blockLength=20, deny=deny, denyOption="trim")[[1]]

  ## suppressPackageStartupMessages(library(BentoBox))
  ## plotGRanges <- function(gr) {
  ##   bb_pageCreate(width = 5, height = 2, xgrid = 0,
  ##                 ygrid = 0, showGuides = TRUE)
  ##   chromend <- 100
  ##   p <- bb_params(chromstart = 0, chromend = chromend,
  ##                  x = 0.5, width = 4, height = 1.5,
  ##                  at = seq(0, chromend, 10),
  ##                  fill = palette()[2:6])
  ##   prngs <- bb_plotRanges(data = gr, params = p,
  ##                          chrom = "chr1", y = 1.5,
  ##                          just = c("left", "bottom"),
  ##                          colorby=colorby("block"))
  ##   bb_annoGenomeLabel(plot = prngs, params = p, y = 1.6, sequence=FALSE)
  ## }

  ## plotGRanges(x)
  ## plotGRanges(b)
})
