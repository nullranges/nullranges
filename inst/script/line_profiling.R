library(profvis)
profvis({
  res <- seg_bootstrap_granges(seg,g1,L_b,within_chrom=FALSE,coarse = TRUE)
})

profvis({
  seg_bootstrap_granges(seg,g1,L_b,within_chrom=TRUE,proportion_length = FALSE,coarse = TRUE)
})

microbenchmark::microbenchmark(
  ns <- seq_len(length(L_s)),
  ns <- sort(unique(state))
)

