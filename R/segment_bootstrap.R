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
#' @export
seg_bootstrap_granges <- function(seg, x, L_b, within_chrom=TRUE,proportion_length = TRUE, coarse = FALSE) {
  chrom_lens <- seqlengths(x)[seqnames(x)@values]
  chroms <- as.character(seqnames(x)@values)
  if (within_chrom) {
    res <- lapply(chroms, function(chr) {
      L_c <- chrom_lens[[chr]]
      seg0 <- seg[seqnames(seg) == chr]
      seg_length <- seg0 %>%
        group_by(state) %>%
        summarise(Ls = sum(width)) # derive each states length
      L_s <- seg_length$Ls
      x0 <- ranges(x[seqnames(x) == chr])
      r_prime <- seg_bootstrap_iranges_fast(ranges(seg0), x0, seg0$state, L_c, L_s, L_b)
      GRanges(seqnames = chr, ranges = r_prime, seqlengths = chrom_lens)
    })
    res <- do.call(c, res)
  } else {
    seg_length <- seg %>%
      group_by(state) %>%
      summarise(Ls = sum(width)) # derive each states length
    L_s <- seg_length$Ls
    # store the gene digit for join_overlap_inner
    x$key.x<-seq_len(length(x)) # or check whether have this column in x
    
    res<-seg_bootstrap_iranges(ranges(seg),x,seg$state,seqnames(seg),L_s,L_b,chrom_lens)
  }
  
  return(res)
}
#' @export
seg_bootstrap_iranges_fast <- function(seg, x, state, L_c, L_s, L_b, proportion_length = TRUE, coarse = FALSE) {
  if (proportion_length) {
    # number of states
    ns <- sort(unique(state)) # some chr may lack some states  
    
    # the block width for each segmentation state is scaled
    # down to the segmentation state size, e.g. if segmentation state
    # is half of the chromosome, then the block width is half of L_b
    L_b0 <- round(L_b * L_s / L_c)
    obj <- lapply(ns, function(m) { # loop over segmentation states
      # the length of the block for this state
      L_bs <- L_b0[which(ns==m)]
      # the segmentation for this state
      seg2 <- seg[state == m]
      # fraction that each range of the segmentation comprises of the whole
      p <- (width(seg2)) / sum(width(seg2))
      # number of blocks within each range of the segmentation
      times <- ceiling(L_c / L_b * p)
      j<-which((end(seg2)-L_bs)<start(seg2))
      width(seg2)[j]=L_bs
      # create random start positions within each segment
      random_start <- unlist(mapply(function(time, x, y) runif(time, x, y), times, start(seg2), end(seg2)-L_bs+1))
      
      # shuffle the blocks
      random_start <- ifelse(length(random_start)==1,round(random_start),sample(random_start))
      # create the random blocks
      random_blocks0 <- IRanges(start = random_start, width = L_bs)
      # the positions of the rearranged blocks in this segmentation state
      start_order <- unlist(mapply(function(x, y) seq(from = x, to = y, by = L_bs), start(seg2), end(seg2)))
      
      rearranged_blocks0 <- IRanges(start = start_order, width = L_bs)
      return(list(random_blocks0, rearranged_blocks0))
    })
    obj <- unlist(obj)
    random_blocks <- do.call(c, obj[seq(1, length(obj), by = 2)])
    rearranged_blocks <- do.call(c, obj[seq(2, length(obj), by = 2)])
    block_shift <- start(rearranged_blocks) - start(random_blocks)
    fo <- findOverlaps(random_blocks, x)
    x_prime <- IRanges::shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
    if (coarse) {
      obj_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_c] # faster
    } else {
      obj_prime <- join_overlap_intersect(x_prime, seg) # accurate
    }
  } else {
    # number of blocks within each range of the segmentation
    times <- ceiling((width(seg)) / L_b)
    # create random start positions within each segment
    random_start <- mapply(function(time, x, y) runif(time, x, y), times, start(seg), end(seg))
    # the positions of the rearranged blocks in this segmentation state
    start_order <- mapply(function(x, y) seq(from = x, to = y, by = L_b), start(seg), end(seg))
    obj <- lapply(1:ns, function(m) {
      # shuffle the blocks
      random_start0 <- sample(unlist(random_start[state == m]))
      # create the random blocks
      random_blocks0 <- IRanges(start = random_start0, width = L_b)
      
      start_order0 <- unlist(start_order[state == m])
      rearranged_blocks0 <- IRanges(start = start_order0, width = L_b)
      return(list(random_blocks0, rearranged_blocks0))
    })
    obj <- unlist(obj)
    random_blocks <- do.call(c, obj[seq(1, length(obj), by = 2)])
    rearranged_blocks <- do.call(c, obj[seq(2, length(obj), by = 2)])
    block_shift <- start(rearranged_blocks) - start(random_blocks)
    
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
seg_bootstrap_iranges <- function(seg, x,state,chrnames,L_s, L_b, chrom_lens) { 
  L_c <- sum(chrom_lens)
# the block width for each segmentation state is scaled
# down to the segmentation state size, e.g. if segmentation state
# is half of the chromosome, then the block width is half of L_b
  L_b0 <- round(L_b * L_s / L_c)
  ns <- length(L_s)
  obj <- lapply(1:ns, function(m) { # loop over segmentation states
    # the length of the block for this state
    L_bs <- L_b0[m]
    # the segmentation for this state
    seg2 <- seg[state == m]
    # fraction that each range of the segmentation comprises of the whole
    p <- (width(seg2)) / sum(width(seg2))
    # number of blocks within each range of the segmentation
    times <- ceiling(L_c / L_b * p)
    seqnames = rep(as.character(chrnames[state==m]),times)
    j<-which((end(seg2)-L_bs)<start(seg2))
    width(seg2)[j]=L_bs
    # create random start positions within each segment
    random_start <- unlist(mapply(function(time, x, y) round(runif(time, x, y)), times, start(seg2), end(seg2)-L_bs+1)) # hit issue negative "end(seg2)-L_bs"
    # create the random blocks
    random_blocks_r <- IRanges(start = random_start, width = L_bs)

    # the positions of the rearranged blocks in this segmentation state
    start_order <- unlist(mapply(function(x, y) seq(from = x, to = y, by = L_bs), start(seg2), end(seg2)))
    rearranged_blocks <- IRanges(start = start_order, width = L_bs)
  
    # shuffle the blocks
    index<-sample(length(random_start))
    block_shift <- start(rearranged_blocks)[index] - start(random_blocks_r)
    return(list(random_blocks_r,block_shift,index,seqnames))
  })
  obj <- do.call(c,obj)
  random_blocks_r <- do.call(c, obj[seq(1, length(obj), by = 4)])
  block_shift <- do.call(c, obj[seq(2, length(obj), by = 4)])
  index <- do.call(c, obj[seq(3, length(obj), by = 4)])
  seqnames <- do.call(c, obj[seq(4, length(obj), by = 4)])
  random_blocks <- GRanges(seqnames = seqnames,ranges = random_blocks_r, seqlengths = chrom_lens) # warning message out-of-bound if not use "end(seg2)-L_bs"
  random_blocks$key.random<-seq_len(length(random_blocks))
  
  fo <- join_overlap_inner(x,random_blocks) # Question: do we want to only select once for each gene?

  x_prime <- IRanges::shift(ranges(x[fo$key.x]), block_shift[fo$key.random])

  # keep <- start(x_prime) >= 1 & end(x_prime) <= chrom_lens[seqnames[index[fo$key.random]]] # fast but discard many genes hit boundary
  gr_prime<-GRanges(seqnames =seqnames[index[fo$key.random]],ranges=x_prime)
  chrom_blocks <- GRanges(seqnames = names(chrom_lens),ranges =IRanges(start=rep(1,length(chrom_lens)),width = chrom_lens)) #Question, do we want that as input
  res<-join_overlap_intersect(gr_prime,chrom_blocks)
  return(res)
}

