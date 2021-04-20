library(profvis)
profvis({
  res <- seg_bootstrap_granges(seg,gr,L_b,within_chrom=TRUE,proportion_length = FALSE,coarse = TRUE)
})
profvis({
  res <- seg_bootstrap_granges(seg,gr,L_b,within_chrom=TRUE,proportion_length = TRUE,coarse = TRUE)
})
profvis({
  res_wochr <- seg_bootstrap_granges(seg,gr,L_b,within_chrom=FALSE,proportion_length = FALSE,coarse = TRUE)
})
profvis({
  res_wochr <- seg_bootstrap_granges(seg,gr,L_b,within_chrom=FALSE,proportion_length = TRUE,coarse = TRUE)
})

microbenchmark::microbenchmark(
  seqname <- rep(chroms,times),
  seqname <- unlist(lapply(seq_len(length(chroms)),function(j) rep(chroms[j],times[j])),use.names = FALSE)
)

seg.length<-function(chrnames,times,ns, state){
  seqnames = lapply(seq_len(length(as.character(chrnames))),function(t) rep(as.character(chrnames)[t],times[t]))
  obj <- lapply(ns, function(m) {
    seqnames <- unlist(seqnames[state == m],use.names = FALSE)
    return(list(seqnames))
  })
}

seg.length2<-function(chrnames,times,ns, state){
  seqnames3 <- rep(as.character(chrnames),times)
  names(seqnames3) <- rep(state,times)
  obj <- lapply(ns, function(m) {
    seqnames <- unname(seqnames3[which(names(seqnames3)== m)])
  return(list(seqnames))
  })
}
microbenchmark::microbenchmark(
  seg.length(chrnames,times,ns, state),
  seg.length2(chrnames,times,ns, state)
)
