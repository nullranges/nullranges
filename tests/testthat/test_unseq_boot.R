library(nullranges)
test_that("unsegmented bootstrap works", {

  library(GenomicRanges)
  seq_nms <- rep(c("chr1","chr2","chr3"),c(4,5,2))
  gr <- GRanges(seqnames=seq_nms,
                IRanges(start=c(1,101,121,201,
                                101,201,216,231,401,
                                1,101),
                        width=c(20, 5, 5, 30,
                                20, 5, 5, 5, 30,
                                80, 40)),
                seqlengths=c(chr1=300,chr2=450,chr3=200),
                chrom=as.integer(factor(seq_nms)))

  blockLength <- 100
  gr_prime <- bootRanges(gr, blockLength=blockLength, type="permute", withinChrom=TRUE)
  gr_prime <- bootRanges(gr, blockLength=blockLength, type="bootstrap", withinChrom=TRUE)
  gr_prime <- bootRanges(gr, blockLength=blockLength, type="permute", withinChrom=FALSE)
  gr_prime <- bootRanges(gr, blockLength=blockLength, type="bootstrap", withinChrom=FALSE)

  # empty test
  expect_true(length(gr_prime) > 0)
  
})
