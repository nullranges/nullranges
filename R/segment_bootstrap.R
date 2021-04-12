#' Segmentated bootstrap GRanges
#'
#' @param seg segmentation GRanges a column("state") indicate segmentation states
#' @param x the input GRanges
#' @param L_c the length of the chromosome
#' @param L_s a vector of the length of each segmentation region
#' @param L_b the length of the block
#' @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
#' @param coarse simply check the chromosome ranges
#' 
#' @importFrom dplyr summarise group_by %>%
#'
#' @export
seg_bootstrap_granges <- function(seg, x, L_b, within_chrom=TRUE,proportion_length = TRUE, coarse = FALSE) {
  chrom_lens <- seqlengths(x)[seqnames(x)@values]
  chroms <- as.character(seqnames(x)@values)
  if (within_chrom) {
    obj <- lapply(chroms, function(chr) {
      L_c <- chrom_lens[[chr]]
      seg0 <- seg[seqnames(seg) == chr]
      seg_length <- seg0 %>%
        group_by(state) %>%
        summarise(Ls = sum(width)) # derive each states length
      L_s <- seg_length$Ls
      x0 <- ranges(x[seqnames(x) == chr])
      r_prime <- seg_bootstrap_iranges_fast(ranges(seg0), x0, seg0$state, L_c, L_s, L_b,proportion_length,coarse)
      r_prime
    })
    times <- lengths(obj)
    seqname <- unlist(lapply(seq_len(length(chroms)),function(j) rep(chroms[j],times[j])),use.names = FALSE)
    obj <- do.call(c, obj)
    res <- GRanges(seqnames = seqname, ranges = obj, seqlengths = chrom_lens)
  } else {
    seg_length <- seg %>%
      group_by(state) %>%
      summarise(Ls = sum(width)) # derive each states length
    L_s <- seg_length$Ls
    # store the gene digit for join_overlap_inner
    x$key.x<-seq_len(length(x)) # or check whether have this column in x
    res<-seg_bootstrap_iranges(ranges(seg),x,seg$state,seqnames(seg),L_s,L_b,chrom_lens,proportion_length,coarse)
  }
  return(res)
}

#' @export
seg_bootstrap_iranges_fast <- function(seg, x, state, L_c, L_s, L_b, proportion_length = TRUE, coarse = FALSE) {
  if (proportion_length) {
    ## number of states
    ns <- sort(unique(state)) # some chr may lack some states  
    
    ## the block width for each segmentation state is scaled
    ## down to the segmentation state size, e.g. if segmentation state
    ## is half of the chromosome, then the block width is half of L_b
    L_b0 <- round(L_b * L_s / L_c)
    obj <- lapply(ns, function(m) { # loop over segmentation states
      ## the length of the block for this state
      L_bs <- L_b0[which(ns==m)]
      ## the segmentation for this state
      seg2 <- seg[state == m]
      ## fraction that each range of the segmentation comprises of the whole
      p <- (width(seg2)) / sum(width(seg2))
      ## number of blocks within each range of the segmentation
      times <- ceiling(L_c / L_b * p)
      j<-which((end(seg2)-L_bs)<start(seg2))
      width(seg2)[j]=L_bs
      ## create random start positions within each segment
      random_start <- unlist(lapply(seq_len(length(times)), function(j){
        runif(times[j],start(seg2)[j],end(seg2)[j]-L_bs+1)
      }),use.names = FALSE)
      ## shuffle the blocks
      random_start <- round(sample(random_start,length(random_start)))
      ## the positions of the rearranged blocks in this segmentation state
      start_order <- unlist(lapply(seq_len(length(times)), function(j){
        seq(from = start(seg2)[j], to = end(seg2)[j], by = L_bs)
      }),use.names = FALSE)
      return(list(random_start, start_order))
    })
    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    block_shift <- random_blocks_start - rearranged_blocks_start
    width <- lengths(random_blocks_start0)
    random_blocks <- IRanges(start = random_blocks_start, 
                             width = unlist(lapply(seq_len(length(L_b0)),function(j) rep(L_b0[j],width[j])),use.names = FALSE))
    fo <- findOverlaps(random_blocks, x)
    x_prime <- IRanges::shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
    if (coarse) {
      obj_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_c] # faster
    } else {
      obj_prime <- join_overlap_intersect(x_prime, seg) # accurate
    }
  } else {
    ## number of blocks within each range of the segmentation
    times <- ceiling((width(seg)) / L_b)
    ## create random start positions within each segment
    random_start <- unlist(lapply(seq_len(length(times)), function(j){
      runif(times[j],start(seg)[j],max(start(seg)[j],end(seg)[j]-L_b+1))
    }),use.names = FALSE)
    ## the positions of the rearranged blocks in this segmentation state
    start_order <- unlist(lapply(seq_len(length(times)), function(j){
      seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    }),use.names = FALSE)
    obj <- lapply(ns, function(m) {
      random_start0 <- sample(unlist(random_start[state == m],use.names = FALSE))
      start_order0 <- unlist(start_order[state == m],use.names = FALSE)
      return(list(random_start0, start_order0))
    })

    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    block_shift <- random_blocks_start - rearranged_blocks_start
    random_blocks <- IRanges(start = random_blocks_start, 
                             width = rep(L_b,length(rearranged_blocks_start)))
    fo <- findOverlaps(random_blocks, x)
    x_prime <- IRanges::shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
    if (coarse) {
      obj_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_c] # faster
    } else {
      obj_prime <- join_overlap_intersect(x_prime, seg) # accurate
    }
  }
  
  return(obj_prime)
}


#' @export
seg_bootstrap_iranges <- function(seg, x,state,chrnames,L_s, L_b, chrom_lens, 
                                  proportion_length = TRUE, coarse = FALSE) { 
  if (proportion_length) {
  L_c <- sum(chrom_lens)
 ## the block width for each segmentation state is scaled
 ## down to the segmentation state size, e.g. if segmentation state
 ## is half of the chromosome, then the block width is half of L_b
  L_b0 <- round(L_b * L_s / L_c)
  ns <- seq_len(length(L_s))
  obj <- lapply(ns, function(m) { # loop over segmentation states
    ## the length of the block for this state
    L_bs <- L_b0[which(ns==m)]
    ## the segmentation for this state
    seg2 <- seg[state == m]
    ## fraction that each range of the segmentation comprises of the whole
    p <- (width(seg2)) / sum(width(seg2))
    # number of blocks within each range of the segmentation
    times <- ceiling(L_c / L_b * p)
    seqnames = rep(as.character(chrnames[state==m]),times)
    j<-which((end(seg2)-L_bs)<start(seg2))
    width(seg2)[j]=L_bs
    # create random start positions within each segment
    random_start <- round(unlist(lapply(seq_len(length(times)), function(j){
      runif(times[j],start(seg2)[j],end(seg2)[j]-L_bs+1)
    }),use.names = FALSE))
    ## shuffle the blocks
    # index<-sample(length(random_start))
    # # create the random blocks
    # random_blocks_r <- IRanges(start = random_start, width = L_bs)

    # the positions of the rearranged blocks in this segmentation state
    start_order <- unlist(lapply(seq_len(length(times)), function(j){
      seq(from = start(seg2)[j], to = end(seg2)[j], by = L_bs)
    }),use.names = FALSE)
    # rearranged_blocks <- IRanges(start = start_order, width = L_bs)
  
    # shuffle the blocks
    index<-sample(length(random_start))
    block_shift <- random_start[index] - start_order
    # block_shift <- start(rearranged_blocks)[index] - start(random_blocks_r)
    # return(list(random_blocks_r,block_shift,index,seqnames))
    return(list(random_start,block_shift,seqnames,index))
  })
  random_blocks_start0 <- lapply(obj, `[[`, 1)
  random_blocks_start <- do.call(c, random_blocks_start0)
  block_shift <- do.call(c, lapply(obj, `[[`, 2))
  seqnames <- do.call(c, lapply(obj, `[[`, 3))
  index <- do.call(c, lapply(obj, `[[`, 4))
  ## deal with different segmentation state index 
  width <- lengths(random_blocks_start0)
  add <- c(0,cumsum(width)[-3])
  add <- unlist(lapply(ns, function(j) rep(add[j],width[j])), use.names = FALSE)
  index <- index + add
  # block_shift <- random_blocks_start - rearranged_blocks_start
  
  
  # obj <- do.call(c,obj)
  # random_blocks_r <- do.call(c, obj[seq(1, length(obj), by = 4)])
  # block_shift <- do.call(c, obj[seq(2, length(obj), by = 4)])
  # index <- do.call(c, obj[seq(3, length(obj), by = 4)])
  # seqnames <- do.call(c, obj[seq(4, length(obj), by = 4)])
  random_blocks_r <- IRanges(start = random_blocks_start[index], 
                           width = unlist(lapply(seq_len(length(L_b0)),function(j) rep(L_b0[j],width[j])),use.names = FALSE))
  
  random_blocks <- GRanges(seqnames = seqnames[index],ranges = random_blocks_r, seqlengths = chrom_lens) # warning message out-of-bound if not use "end(seg2)-L_bs"
  random_blocks$key.random<-seq_len(length(random_blocks))
  
  fo <- join_overlap_inner(x,random_blocks) # Question: do we want to only select once for each gene?

  x_prime <- IRanges::shift(ranges(x[fo$key.x]), block_shift[fo$key.random])
  if (coarse) {
    obj_prime <- start(x_prime) >= 1 & end(x_prime) <= chrom_lens[seqnames[fo$key.random]] # faster
    res <- GRanges(seqnames = seqnames[fo$key.random][obj_prime],ranges=x_prime[obj_prime],seqlengths = chrom_lens)
  } else {
    suppressWarnings(gr_prime<-GRanges(seqnames = seqnames[fo$key.random],ranges=x_prime,seqlengths = chrom_lens))
    res <-trim(gr_prime)
  }
  }else {
    ## number of blocks within each range of the segmentation
    times <- ceiling((width(seg)) / L_b)
    seqnames = rep(as.character(chrnames),times)
    ## create random start positions within each segment
    random_start <- unlist(lapply(seq_len(length(times)), function(j){
      runif(times[j],start(seg)[j],max(start(seg)[j],end(seg)[j]-L_b+1))
    }),use.names = FALSE)
    ## the positions of the rearranged blocks in this segmentation state
    start_order <- unlist(lapply(seq_len(length(times)), function(j){
      seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    }),use.names = FALSE)
    obj <- lapply(ns, function(m) {
      random_start0 <- sample(unlist(random_start[state == m],use.names = FALSE))
      start_order0 <- unlist(start_order[state == m],use.names = FALSE)
      return(list(random_start0, start_order0))
    })
    
    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    block_shift <- random_blocks_start - rearranged_blocks_start
    random_blocks <- IRanges(start = random_blocks_start, 
                             width = rep(L_b,length(rearranged_blocks_start)))
    fo <- findOverlaps(random_blocks, x)
    x_prime <- IRanges::shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
    if (coarse) {
      obj_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_c] # faster
    } else {
      obj_prime <- join_overlap_intersect(x_prime, seg) # accurate
    }
  }
  return(res)
}

