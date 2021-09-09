library(nullranges)
test_that("trim excludeOption works as expected", {

  x <- GRanges("chr1", IRanges(1 + 2 * 0:45, width=10),
               strand=rep(c("+","-"),length=46),
               seqlengths=c(chr1=100))
  exclude <- GRanges("chr1", IRanges(36,45),
                     seqlengths=c(chr1=100))
  x <- sort(x)
  x$block <- 1
  b <- bootRanges(x, blockLength=20, exclude=exclude, excludeOption="trim")[[1]]

})
