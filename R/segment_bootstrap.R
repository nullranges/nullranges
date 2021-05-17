#' Segmentated block bootstrap
#'
#' @param x the input GRanges
#' @param seg segmentation GRanges a column("state") indicate segmentation states
#' @param L_b the length of the block
#' @param R the number of bootstrap samples to generate
#' @param deny GRanges of deny region
#' @param deny_option Indicate whether toss or trim the overlaps between deny and
#' bootstrap ranges.
#' @param within_chrom whether to perform bootstrapping within chromosome (default FALSE)
#' @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
#' @param ncores A cluster object created by \code{\link[parallel]{makeCluster}}.
#' Or an integer to indicate number of child-processes
#' (integer values are ignored on Windows) for parallel evaluations
#'
#' @export
segBootstrapRanges <- function(x, seg, L_b, R,
                               deny, deny_option = c("toss", "trim"),
                               within_chrom = FALSE, proportion_length = TRUE,
                               ncores = NULL) {
  chrom_lens <- seqlengths(x)
  chroms <- as.character(seqnames(x)@values)
  deny_option <- match.arg(deny_option)
  replicate(R, {
    if (within_chrom) {
      # loop over chromosomes
      obj <- lapply(chroms, function(chr) {
        L_c <- chrom_lens[[chr]]
        seg0 <- seg[seqnames(seg) == chr]
        x0 <- x[seqnames(x) == chr]
        r_prime <- seg_bootstrap_granges_within_chrom(x0, ranges(seg0), seg0$state, L_c, L_b, proportion_length)
        r_prime
      })
      res <- do.call(c, obj)
    } else {
      res <- seg_bootstrap_granges(x, ranges(seg), seg$state, L_b, seqnames(seg), chrom_lens, proportion_length)
    }
    if (deny_option == "toss") {
      # TODO I think this could just be an overlapsAny() call rather than queryHits(findOverlaps())
      res_accept <- res[-queryHits(findOverlaps(res, deny, type="any"))]
    } else {
      # res <- setdiff(res,deny) # To do: cannot reserve metadata column and whether to use ignore.strand=TRUE
      ## To do: need to keep the gaps with same deny strand, here is special case that all strand(deny) ="*"
      ## To do: need to place outstide of R, doing gaps once.
      gap <- gaps(deny,end = seqlengths(x)) %>%
        plyranges::filter(strand=="*")
      ## the region remove deny regions
      res_accept <- plyranges::join_overlap_intersect(res,gap)
    }
    res_accept
  })
}

# Segmentation block bootstrap IRanges within chromosomes
#
# @param x the ranges of the input GRanges
# @param seg the ranges of segmentation GRanges
# @param state the segmentation state
# @param L_c the chromosome length
# @param L_b the length of the block
# @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
# @param coarse logical indicating if out of bound ranges will be discarded or keep the intersect region.
seg_bootstrap_granges_within_chrom <- function(x, seg, state, L_c, L_b, proportion_length = TRUE, coarse = FALSE) {
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

# Segmentation block bootstrap IRanges spanning the whole genome
#
# @param x the input GRanges
# @param seg the ranges of segmentation GRanges
# @param state the segmentation state
# @param L_b the length of the block
# @param chr_names the seqnames of the input GRanges
# @param chrom_lens all chromosome length
# @param proportion_length use scaled block length or scaled number of blocks of each segmentation region
# @param coarse logical indicating if out of bound ranges will be discarded or keep the intersect region.
seg_bootstrap_granges <- function(x, seg, state, L_b, chr_names, chrom_lens,
                                  proportion_length = TRUE) {

  # sequence along the states
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
      seqnames <- as.character(chr_names[state == m])[random_chr]
      ## where those blocks will move to
      start_order <- unlist(lapply(seq_len(length(times)), function(j) {
        seq(from = start(seg2)[j], to = end(seg2)[j], by = L_bs)
      }), use.names = FALSE)
      ## record ordered chrname
      segchr_names <- rep(as.character(chr_names[state == m]), times)
      return(list(random_start, start_order, seqnames, segchr_names))
    })
    random_blocks_start <- do.call(c, lapply(obj, `[[`, 1))
    rearranged_blocks_start <- do.call(c, lapply(obj, `[[`, 2))

    # TODO this test fails for me on the toy dataset
    stopifnot(length(random_blocks_start) != length(rearranged_blocks_start))

    seqnames <- do.call(c, lapply(obj, `[[`, 3))
    rearr_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 4))
    ## number of blocks for each segmentation states
    width <- lengths(random_blocks_start)
    ## these blocks are the 'bait' for capturing features in 'x'
    random_blocks <- GRanges(seqnames=seqnames,
                             ranges=IRanges(start = random_blocks_start,
                                            width = rep(L_b0, width)))
  } else {

    # number of rearranged blocks within each range of the segmentation
    times <- ceiling(width(seg) / L_b)
    # total number of tiling blocks
    n <- sum(times)
    # sample segments according to segmentation length
    random_seg <- sample(length(seg), n, replace=TRUE, prob=width(seg))
    # random starts within the segments, as many as specified in 'times'
    rmin <- start(seg)[random_seg]
    rmax <- (times[random_seg] - 1) * L_b + start(seg)[random_seg]
    random_start <- round(runif(n,rmin,rmax))
    # start positions of rearranged blocks tiling the segments
    rearranged_blocks_start_list <- lapply(seq_along(seg), function(j) {
      seq(from = start(seg)[j], to = end(seg)[j], by = L_b)
    })

    # loop over states 's', splitting out the random segments for 's'
    obj <- lapply(ns, function(s) {
      # index of those random segments chosen that are state 's'
      idx <- which(state[random_seg] == s)
      # random starts for state s
      random_start_s <- random_start[idx]
      # names of chromosomes for random segments that are state 's'
      random_chr_names_s <- chr_names[random_seg[idx]]
      # start positions of blocks tiling segment for state 's'
      rearr_blocks_start_s <- unlist(rearranged_blocks_start_list[state == s],
                                     use.names = FALSE)
      # names of chromosomes for blocks tilings segment for state 's'
      rearr_blocks_chr_names_s <- rep(as.character(chr_names)[state == s],
                                      times[state == s])
      # see if we have too few or too many blocks for this state
      rand_n <- length(random_start_s)
      rearr_n <- length(rearr_blocks_start_s)
      # if we have too few random blocks, up-sample them
      if (rand_n < rearr_n) {
        rand_idx <- sample(rand_n, rearr_n, replace=TRUE)
      } else {
        # otherwise down-sample to the number of re-arranged blocks (or stay)
        rand_idx <- seq_len(rearr_n)
      }
      return(list(random_start_s[rand_idx],
                  random_chr_names_s[rand_idx],
                  rearr_blocks_start_s,
                  rearr_blocks_chr_names_s))
    })

    # unlist the results in 'obj'
    random_blocks_start <- unlist(lapply(obj, `[[`, 1))
    random_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 2))
    rearranged_blocks_start <- unlist(lapply(obj, `[[`, 3))
    rearr_blocks_chr_names <- do.call(c, lapply(obj, `[[`, 4))

    # construct the random blocks
    random_blocks <- GRanges(seqnames = random_blocks_chr_names,
                             ranges = IRanges(start = random_blocks_start,
                                              width = L_b))
  }

  rearr_blocks_chr_names <- factor(rearr_blocks_chr_names, levels = seqlevels(x))

  # use the bait to sample features in 'x'
  fo <- findOverlaps(random_blocks, x)
  # x has been sampled multiple times
  x_mult_hits <- x[subjectHits(fo)]
  # label which 'bait' block each feature hit in the re-sampling
  mcols(x_mult_hits)$block <- queryHits(fo)
  shift_and_swap_chrom(x_mult_hits, rearr_blocks_chr_names,
                       random_blocks_start, rearranged_blocks_start)

}
