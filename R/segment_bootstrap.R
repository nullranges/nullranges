#' Segmentated bootstrap GRanges
#'
#' @param seg segmentation GRanges a column("state") indicate segmentation states
#' @param x the input GRanges
#' @param L_b the length of the block
#' @param deny GRanges of deny region
#' @param R the number time of bootstrap
#' @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
#' @param ncores A cluster object created by \code{\link[parallel]{makeCluster}}.
#' Or an integer to indicate number of child-processes
#' (integer values are ignored on Windows) for parallel evaluations
#'
#' @importFrom pbapply pblapply
#'
#' @export
segBootstrapRanges <- function(seg, x, L_b, deny, R, within_chrom = FALSE,
                                  proportion_length = TRUE, ncores = NULL) {
  chrom_lens <- seqlengths(x)
  chroms <- as.character(seqnames(x)@values)
  ans <- pblapply(seq_len(R), function(i) {
    if (within_chrom) {
      obj <- lapply(chroms, function(chr) {
        L_c <- chrom_lens[[chr]]
        seg0 <- seg[seqnames(seg) == chr]
        x0 <- x[seqnames(x) == chr]
        r_prime <- seg_bootstrap_granges_within_chrom(ranges(seg0), x0, seg0$state, L_c, L_b, proportion_length)
        r_prime
      })
      res <- do.call(c, obj)
    } else {
      res <- seg_bootstrap_granges(ranges(seg), x, seg$state, L_b, seqnames(seg), chrom_lens, proportion_length)
    }
    res <- setdiff(res,deny) # To do: whether to use ignore.strand=TRUE ? 
    return(res)
  }, cl = ncores)
}

# segmentation Block bootstrap IRanges within chromosome
#
# @param seg the ranges of segmentation GRanges
# @param x the ranges of the input GRanges
# @param state the segmentation state
# @param L_c the chromosome length
# @param L_b the length of the block
# @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
# @param coarse logical indicating if out of bound ranges will be discarded or keep the intersect region.
seg_bootstrap_granges_within_chrom <- function(seg, x, state, L_c, L_b, proportion_length = TRUE, coarse = FALSE) {
  ## number of states
  ns <- sort(unique(state)) # some chr may lack some states
  if (proportion_length) {
    ## Derive segmentation state length
    L_s <- unlist(lapply(ns, function(j) sum(width(seg)[which(state == j)])), use.names = FALSE)

    ## the block width for each segmentation state is scaled
    ## down to the segmentation state size, e.g. if segmentation state
    ## is half of the chromosome, then the block width is half of L_b
    L_b0 <- round(L_b * L_s / L_c)
    obj <- lapply(ns, function(m) { # loop over segmentation states
      ## the length of the block for this state
      L_bs <- L_b0[which(ns == m)]
      ## the segmentation for this state
      seg2 <- seg[state == m]
      ## fraction that each range of the segmentation comprises of the whole
      p <- (width(seg2)) / sum(width(seg2))
      ## number of blocks within each range of the segmentation
      times <- ceiling(L_c / L_b * p)
      j <- which((end(seg2) - L_bs) < start(seg2))
      width(seg2)[j] <- L_bs
      ## create random start positions within each segment
      random_start <- unlist(lapply(seq_len(length(times)), function(j) {
        runif(times[j], start(seg2)[j], end(seg2)[j] - L_bs + 1)
      }), use.names = FALSE)
      ## shuffle the blocks
      random_start <- round(sample(random_start, length(random_start)))
      ## the positions of the rearranged blocks in this segmentation state
      start_order <- unlist(lapply(seq_len(length(times)), function(j) {
        seq(from = start(seg2)[j], to = end(seg2)[j], by = L_bs)
      }), use.names = FALSE)
      return(list(random_start, start_order))
    })
    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    block_shift <- random_blocks_start - rearranged_blocks_start
    width <- lengths(random_blocks_start0)
    random_blocks <- IRanges(
      start = random_blocks_start,
      width = rep(L_b0, width)
    )
  } else {
    ## number of blocks within each range of the segmentation
    times <- ceiling((width(seg)) / L_b)
    n <- sum(times)
    index <- seq_len(length(times))
    # sample according to segmentation length
    random_chr <- sample(index, n, replace=TRUE, prob=width(seg))
    # random starts within the tiled regions
    random_start <- round(runif(n, start(seg)[random_chr], (times[random_chr] -1 ) * L_b + start(seg)[random_chr]))
    ## the positions of the rearranged blocks in this segmentation state
    start_order <- lapply(seq_len(length(times)), function(j) {
      seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    })
    obj <- lapply(ns, function(m) {
      poi <- which(state == m)
      index <- which(random_chr %in% poi)
      random_start0 <- random_start[index]
      start_order0 <- unlist(start_order[state == m], use.names = FALSE)
      return(list(random_start0, start_order0))
    })
    random_blocks_start <- do.call(c,lapply(obj, `[[`, 1))
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    block_shift <- random_blocks_start - rearranged_blocks_start
    random_blocks <- IRanges(start = random_blocks_start,width = L_b)
  }
  fo <- findOverlaps(random_blocks, ranges(x))
  # shift the ranges in those bait blocks
  suppressWarnings({x_prime <- shift(x[subjectHits(fo)], block_shift[queryHits(fo)])})
  x_prime <- trim(x_prime)
  return(x_prime)
}

# segmentation Block bootstrap IRanges span the whole genome.
#
# @param seg the ranges of segmentation GRanges
# @param x the input GRanges with a column ("key.x") indicate the index of gene
# @param state the segmentation state
# @param L_b the length of the block
# @param chrnames the seqnames of the input GRanges
# @param chrom_lens all chromosome length
# @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
# @param coarse logical indicating if out of bound ranges will be discarded or keep the intersect region.
seg_bootstrap_granges <- function(seg, x, state, L_b, chrnames, chrom_lens,
                                   proportion_length = TRUE) {
  ns <- seq_len(max(state))
  if (proportion_length) {
    L_c <- sum(chrom_lens)
    ## Derive segmentation state length
    L_s <- unlist(lapply(seq_len(max(state)), function(j) sum(width(seg)[which(state == j)])), use.names = FALSE)
    ## the block width for each segmentation state is scaled
    ## down to the segmentation state size, e.g. if segmentation state
    ## is half of the chromosome, then the block width is half of L_b
    L_b0 <- round(L_b * L_s / L_c)
    obj <- lapply(ns, function(m) { # loop over segmentation states
      ## the length of the block for this state
      L_bs <- L_b0[which(ns == m)]
      ## the segmentation for this state
      seg2 <- seg[state == m]
      ## fraction that each range of the segmentation comprises of the whole
      p <- (width(seg2)) / sum(width(seg2))
      ## number of blocks within each range of the segmentation
      times <- ceiling(L_c / L_b * p)
      ## total number of tiling blocks
      n <- sum(times)
      index <- seq_len(length(times))
      ## sample according to segmentation length
      random_chr <- sample(index, n, replace=TRUE, prob=width(seg2))
      ## random starts within the tiled regions
      random_start <- round(runif(n, start(seg2)[random_chr], (times[random_chr] -1 ) * L_bs + start(seg2)[random_chr]))
      ## record sampled chr name
      seqnames <- as.character(chrnames[state == m])[random_chr]
      ## where those blocks will move to
      start_order <- unlist(lapply(seq_len(length(times)), function(j) {
        seq(from = start(seg2)[j], to = end(seg2)[j], by = L_bs)
      }), use.names = FALSE)
      ## record ordered chrname
      segchrnames <- rep(as.character(chrnames[state == m]), times)
      return(list(random_start, start_order, seqnames, segchrnames))
    })
    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    seqnames <- do.call(c, lapply(obj, `[[`, 3))
    segchrnames <- do.call(c, lapply(obj, `[[`, 4))
    ## number of blocks for each segmentation states
    width <- lengths(random_blocks_start0)
    ## these blocks are the 'bait' for capturing features in 'x'
    random_blocks <- GRanges(seqnames=seqnames,
                             ranges=IRanges(start = random_blocks_start, width = rep(L_b0, width)))
  } else {
    ## number of blocks within each range of the segmentation
    times <- ceiling((width(seg)) / L_b)
    # total number of tiling blocks
    n <- sum(times)
    index <- seq_len(length(times))
    # sample according to segmentation length
    random_chr <- sample(index, n, replace=TRUE, prob=width(seg))
    # random starts within the tiled regions
    random_start <- round(runif(n, start(seg)[random_chr], (times[random_chr] -1 ) * L_b + start(seg)[random_chr]))
    # where those blocks will move to
    start_order <- lapply(seq_len(length(times)), function(j) {
      seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    })
    ## record ordered chrname
    segchrnames <- rep(as.character(chrnames), times)
    obj <- lapply(ns, function(m) {
      poi <- which(state == m)
      index <- which(random_chr %in% poi)
      random_start0 <- random_start[index]
      start_order0 <- unlist(start_order[state == m], use.names = FALSE)
      ## record sampled chr name
      seqnames <- chrnames[random_chr[index]]
      return(list(random_start0, start_order0, seqnames))
    })
    random_blocks_start <- do.call(c, lapply(obj, `[[`, 1))
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    ## deal with different segmentation state index
    seqnames <- do.call(c, lapply(obj, `[[`, 3))
    random_blocks <- GRanges(seqnames=seqnames,
                             ranges=IRanges(start = random_blocks_start, width = L_b))
  }
  segchrnames <- factor(segchrnames,levels = seqlevels(x))
  ## use the bait to sample features in 'x'
  fo <- findOverlaps(random_blocks, x)
  ## x has been sampled multiple times
  x_mult_hits <- x[subjectHits(fo)]
  ## label which 'bait' block each feature hit in the re-sampling
  mcols(x_mult_hits)$block <- queryHits(fo)
  res <- shift_and_swap_chrom(x_mult_hits, segchrnames,random_blocks_start, rearranged_blocks_start)
  return(res)
}
