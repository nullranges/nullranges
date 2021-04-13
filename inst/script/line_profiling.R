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
  seg.length(seg,times,ns,state),
  seg.length2(seg,times,ns,state)
)

seg.length<-function(seg,times,ns, state){
  start <- lapply(seq_len(length(times)), function(j){
    random_start <- runif(times[j],start(seg)[j],max(start(seg)[j],end(seg)[j]-L_b+1))
    start_order <- seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    return(list(random_start,start_order))
  })
  obj <- lapply(ns, function(m) {
    poi <- start[state == m]
    random_start0 <- sample(unlist(lapply(poi, `[[`, 1),use.names = FALSE))
    start_order0 <- unlist(lapply(poi, `[[`, 2),use.names = FALSE)
    return(list(random_start0, start_order0))
  })
}

seg.length2<-function(seg,times,ns,state){
  random_start <- lapply(seq_len(length(times)), function(j){
    runif(times[j],start(seg)[j],max(start(seg)[j],end(seg)[j]-L_b+1))
  })
  ## the positions of the rearranged blocks in this segmentation state
  start_order <- lapply(seq_len(length(times)), function(j){
    seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
  })
  obj <- lapply(ns, function(m) {
    poi <- state == m
    random_start0 <- sample(unlist(random_start[poi],use.names = FALSE))
    start_order0 <- unlist(start_order[poi],use.names = FALSE)
    return(list(random_start0, start_order0))
  })
}
microbenchmark::microbenchmark(
  seq_len(length(c(1,2,3))),
  sort(unique(c(1,2,3)))
)
