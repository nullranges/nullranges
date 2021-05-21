context("unsegmented bootstrap")
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

  L_b <- 100
  gr_prime <- bootRanges(gr, L_b=L_b, type="permute", within_chrom=TRUE)
  gr_prime <- bootRanges(gr, L_b=L_b, type="bootstrap", within_chrom=TRUE)
  gr_prime <- bootRanges(gr, L_b=L_b, type="permute", within_chrom=FALSE)
  gr_prime <- bootRanges(gr, L_b=L_b, type="bootstrap", within_chrom=FALSE)

  # empty test
  expect_true(length(gr_prime) > 0)
  
})
