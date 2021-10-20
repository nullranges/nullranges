library(nullranges)
test_that("unseg and seg bootstrap works", {

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
  expect_true(all(table(seqnames(gr)) == table(seqnames(gr_prime[[1]]))))
    
  gr_prime <- bootRanges(gr, blockLength, type="bootstrap", withinChrom=TRUE)
  gr_prime <- bootRanges(gr, blockLength, type="permute", withinChrom=FALSE)
  gr_prime <- bootRanges(gr, blockLength, type="bootstrap", withinChrom=FALSE)

  seg <- GRanges(c("chr1","chr1","chr2","chr2","chr3"),
                 IRanges(c(1,151,1,226,1),c(150,300,225,450,200)),
                 state=c(1,2,1,2,1))

  fo <- findOverlaps(gr, seg)
  gr$seg_state <- seg$state[subjectHits(fo)]
  gr_prime <- bootRanges(gr, blockLength, seg=seg)
  sort(gr_prime)

})
