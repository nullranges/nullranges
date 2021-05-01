#' Segmentated bootstrap GRanges
#'
#' @param seg segmentation GRanges a column("state") indicate segmentation states
#' @param x the input GRanges
#' @param L_c the length of the chromosome
#' @param L_s a vector of the length of each segmentation region
#' @param L_b the length of the block
#' @param R the number time of bootstrap
#' @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
#' @param coarse logical indicating if out of bound ranges will be discarded or keep the intersect region.
#' @param ncores A cluster object created by \code{\link[parallel]{makeCluster}}.
#' Or an integer to indicate number of child-processes
#' (integer values are ignored on Windows) for parallel evaluations
#'
#' @importFrom dplyr summarise group_by %>%
#' @importFrom pbapply pblapply
#'
#' @export
seg_bootstrap_granges <- function(seg, x, L_b, R, within_chrom = TRUE,
                                  proportion_length = TRUE, coarse = FALSE, ncores = NULL) {
  chrom_lens <- seqlengths(x)[seqnames(x)@values]
  chroms <- as.character(seqnames(x)@values)
  ans <- pblapply(seq_len(R), function(i) {
    if (within_chrom) {
      obj <- lapply(chroms, function(chr) {
        L_c <- chrom_lens[[chr]]
        seg0 <- seg[seqnames(seg) == chr]
        x0 <- ranges(x[seqnames(x) == chr])
        r_prime <- seg_bootstrap_iranges_wchr(ranges(seg0), x0, seg0$state, L_c, L_b, proportion_length, coarse)
        r_prime
      })
      times <- lengths(obj)
      seqname <- rep(chroms, times)
      obj <- do.call(c, obj)
      res <- GRanges(seqnames = seqname, ranges = obj, seqlengths = chrom_lens)
    } else {
      # store the gene index for join_overlap_inner
      x$key.x <- seq_len(length(x)) # or check whether have this column in x
      res <- seg_bootstrap_iranges(ranges(seg), x, seg$state, L_b, seqnames(seg), chrom_lens, proportion_length, coarse)
    }
    return(res)
  }, cl = ncores)
}

#' segmentation Block bootstrap IRanges within chromosome
#'
#' @param seg the ranges of segmentation GRanges
#' @param x the ranges of the input GRanges
#' @param state the segmentation state
#' @param L_c the chromosome length
#' @param L_b the length of the block
#' @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
#' @param coarse logical indicating if out of bound ranges will be discarded or keep the intersect region.
#'
#' @export
seg_bootstrap_iranges_wchr <- function(seg, x, state, L_c, L_b, proportion_length = TRUE, coarse = FALSE) {
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
    ## create random start positions within each segment
    random_start <- lapply(seq_len(length(times)), function(j) {
      runif(times[j], start(seg)[j], max(start(seg)[j], end(seg)[j] - L_b + 1))
    })
    ## the positions of the rearranged blocks in this segmentation state
    start_order <- lapply(seq_len(length(times)), function(j) {
      seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    })
    obj <- lapply(ns, function(m) {
      random_start0 <- sample(unlist(random_start[state == m], use.names = FALSE))
      start_order0 <- unlist(start_order[state == m], use.names = FALSE)
      return(list(random_start0, start_order0))
    })

    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    block_shift <- random_blocks_start - rearranged_blocks_start
    random_blocks <- IRanges(
      start = random_blocks_start,
      width = rep(L_b, length(rearranged_blocks_start))
    )
  }
  fo <- findOverlaps(random_blocks, x)
  x_prime <- IRanges::shift(x[subjectHits(fo)], block_shift[queryHits(fo)])
  if (coarse) {
    obj_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_c] # faster
  } else {
    obj_prime <- join_overlap_intersect(x_prime, seg) # accurate
  }
  return(obj_prime)
}

#' segmentation Block bootstrap IRanges span the whole genome.
#'
#' @param seg the ranges of segmentation GRanges
#' @param x the input GRanges with a column ("key.x") indicate the index of gene
#' @param state the segmentation state
#' @param L_b the length of the block
#' @param chrnames the seqnames of the input GRanges
#' @param chrom_lens all chromosome length
#' @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
#' @param coarse logical indicating if out of bound ranges will be discarded or keep the intersect region.
#'
#' @importFrom plyranges join_overlap_inner
#'
#' @export
seg_bootstrap_iranges <- function(seg, x, state, L_b, chrnames, chrom_lens,
                                  proportion_length = TRUE, coarse = FALSE) {
  if (proportion_length) {
    L_c <- sum(chrom_lens)
    ## Derive segmentation state length
    L_s <- unlist(lapply(seq_len(max(state)), function(j) sum(width(seg)[which(state == j)])), use.names = FALSE)
    ## the block width for each segmentation state is scaled
    ## down to the segmentation state size, e.g. if segmentation state
    ## is half of the chromosome, then the block width is half of L_b
    L_b0 <- round(L_b * L_s / L_c)
    ns <- seq_len(length(L_s))
    obj <- lapply(ns, function(m) { # loop over segmentation states
      ## the length of the block for this state
      L_bs <- L_b0[which(ns == m)]
      ## the segmentation for this state
      seg2 <- seg[state == m]
      ## fraction that each range of the segmentation comprises of the whole
      p <- (width(seg2)) / sum(width(seg2))
      # number of blocks within each range of the segmentation
      times <- ceiling(L_c / L_b * p)
      seqnames <- rep(as.character(chrnames[state == m]), times)
      j <- which((end(seg2) - L_bs) < start(seg2))
      width(seg2)[j] <- L_bs
      # create random start positions within each segment
      random_start <- round(unlist(lapply(seq_len(length(times)), function(j) {
        runif(times[j], start(seg2)[j], end(seg2)[j] - L_bs + 1)
      }), use.names = FALSE))

      # the positions of the rearranged blocks in this segmentation state
      start_order <- unlist(lapply(seq_len(length(times)), function(j) {
        seq(from = start(seg2)[j], to = end(seg2)[j], by = L_bs)
      }), use.names = FALSE)

      # shuffle the blocks
      index <- sample(length(random_start))

      return(list(random_start, start_order, seqnames, index))
    })
    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    seqnames <- do.call(c, lapply(obj, `[[`, 3))
    index <- do.call(c, lapply(obj, `[[`, 4))
    ## deal with different segmentation state index
    width <- lengths(random_blocks_start0)
    add <- c(0, cumsum(width)[-length(width)])
    add2 <- rep(add, width)
    index <- index + add

    block_shift <- random_blocks_start[index] - rearranged_blocks_start
    random_blocks_r <- IRanges(
      start = random_blocks_start[index],
      width = rep(L_b0, width)
    )
  } else {
    ns <- seq_len(max(state))
    ## number of blocks within each range of the segmentation
    times <- ceiling((width(seg)) / L_b)
    seqnames <- rep(as.character(chrnames), times)
    names(seqnames) <- rep(state, times)
    ## create random start positions within each segment
    random_start <- lapply(seq_len(length(times)), function(j) {
      runif(times[j], start(seg)[j], max(start(seg)[j], end(seg)[j] - L_b + 1))
    })
    ## the positions of the rearranged blocks in this segmentation state
    start_order <- lapply(seq_len(length(times)), function(j) {
      seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    })
    obj <- lapply(ns, function(m) {
      random_start0 <- unlist(random_start[state == m], use.names = FALSE)
      index <- sample(length(random_start0))
      start_order0 <- unlist(start_order[state == m], use.names = FALSE)
      seqnames <- unname(seqnames[which(names(seqnames) == m)])
      return(list(random_start0, start_order0, index, seqnames))
    })

    random_blocks_start0 <- lapply(obj, `[[`, 1)
    random_blocks_start <- do.call(c, random_blocks_start0)
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))
    index <- do.call(c, lapply(obj, `[[`, 3))
    ## deal with different segmentation state index
    width <- lengths(random_blocks_start0)
    add <- c(0, cumsum(width)[-length(width)])
    add <- rep(add, width)
    index <- index + add
    seqnames <- do.call(c, lapply(obj, `[[`, 4))

    block_shift <- random_blocks_start[index] - rearranged_blocks_start

    random_blocks_r <- IRanges(
      start = random_blocks_start[index],
      width = rep(L_b, length(random_blocks_start))
    )
  }
  ## warning message out-of-bound when end(seg)[j]-L_b+1<0
  suppressWarnings({
    random_blocks <- GRanges(
      seqnames = seqnames[index],
      ranges = random_blocks_r, seqlengths = chrom_lens
    )
  })
  suppressWarnings({
    random_blocks$key.random <- seq_len(length(random_blocks))
  })
  # Question1: do we want to only select once for each gene?
  # Question2: Join_overlap_inner or join_overlap_intersect
  fo <- plyranges::join_overlap_inner(x, random_blocks)
  x_prime <- IRanges::shift(ranges(x[fo$key.x]), block_shift[fo$key.random])
  if (coarse) {
    obj_prime <- start(x_prime) >= 1 & end(x_prime) <=
      chrom_lens[seqnames[fo$key.random]] # faster
    res <- GRanges(
      seqnames = seqnames[fo$key.random][obj_prime],
      ranges = x_prime[obj_prime], seqlengths = chrom_lens
    )
  } else {
    suppressWarnings({
      gr_prime <- GRanges(
        seqnames = seqnames[fo$key.random],
        ranges = x_prime, seqlengths = chrom_lens
      )
    })
    res <- trim(gr_prime)
  }
  return(res)
}
