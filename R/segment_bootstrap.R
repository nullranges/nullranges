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
seg_bootstrap_granges <- function(seg, x, L_b, proportion_length = TRUE, coarse = FALSE) {
  chrom_lens <- seqlengths(x)[seqnames(x)@values]
  chroms <- as.character(seqnames(x)@values)

  seg_length <- seg %>%
    group_by(state) %>%
    summarise(Ls = sum(width)) # derive each states length
  L_s <- seg_length$Ls
  state<-seg$state
  ns <- length(L_s)
  res<-seg_bootstrap_iranges(seg,x,state,L_s,L_b,chrom_lens)
  return(res)
}
seg_bootstrap_iranges <- function(seg, x, state,  L_s, L_b, chrom_lens) {  
  L_c <- sum(chrom_lens)
# the block width for each segmentation state is scaled
# down to the segmentation state size, e.g. if segmentation state
# is half of the chromosome, then the block width is half of L_b
  L_b0 <- round(L_b * L_s / L_c)
  # store the gene digit
  x$key.x<-seq_len(length(x))
  obj <- lapply(1:ns, function(m) { # loop over segmentation states
    # the length of the block for this state
    L_bs <- L_b0[m]
    # the segmentation for this state
    seg2 <- seg[state == m]
    # fraction that each range of the segmentation comprises of the whole
    p <- (width(seg2)) / sum(width(seg2))
    # number of blocks within each range of the segmentation
    times <- ceiling(L_c / L_b * p)
    seqnames = rep(as.character(seqnames(seg2)),times)
    # create random start positions within each segment
    random_start <- unlist(mapply(function(time, x, y) runif(time, x, y), times, start(seg2), end(seg2))) # hit issue negative "end(seg2)-L_bs"
    # create the random blocks
    random_blocks_r <- IRanges(start = random_start, width = L_bs)
    # random_blocks0 <- GRanges(seqnames = seqnames,ranges = random_blocks_r, seqlengths = chrom_lens, # warning message out-of-bound if not use "end(seg2)-L_bs"
    #                         key.random = seq_len(length(random_start)))
    random_blocks0 <- GRanges(seqnames = seqnames,ranges = random_blocks_r, seqlengths = chrom_lens) # warning message out-of-bound if not use "end(seg2)-L_bs"
    # the positions of the rearranged blocks in this segmentation state
    start_order <- unlist(mapply(function(x, y) seq(from = x, to = y, by = L_bs), start(seg2), end(seg2)))
    rearranged_blocks0 <- IRanges(start = start_order, width = L_bs)
  
    # shuffle the blocks
    index<-sample(length(random_start))
    block_shift <- start(rearranged_blocks0)[index] - start(random_blocks_r)
    return(list(random_blocks0,block_shift,index,seqnames))
  })
  obj <- do.call(c,obj)
  random_blocks <- do.call(c, obj[seq(1, length(obj), by = 4)])
  random_blocks$key.random<-seq_len(length(random_blocks))
  block_shift <- do.call(c, obj[seq(2, length(obj), by = 4)])
  index <- do.call(c, obj[seq(3, length(obj), by = 4)])
  seqnames <- do.call(c, obj[seq(4, length(obj), by = 4)])
  
  fo <- join_overlap_inner(x,random_blocks)

  x_prime <- IRanges::shift(ranges(x[fo$key.x]), block_shift[fo$key.random])

  keep <- start(x_prime) >= 1 & end(x_prime) <= chrom_lens[seqnames[index[fo$key.random]]] # warning message out-of-bound
  res<-GRanges(seqnames =seqnames[index[fo$key.random]][keep],ranges=x_prime[keep],seqlengths = chrom_lens)
  
  return(res)
}

