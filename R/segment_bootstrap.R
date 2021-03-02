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
seg_bootstrap_granges <- function(seg, x, within_chrom = TRUE, L_b, proportion_length = TRUE, coarse = FALSE) {
  # type <- match.arg(type)
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
      r_prime <- seg_bootstrap_iranges(ranges(seg0), x0, seg0$state, L_c, L_s, L_b)
      GRanges(seqnames = chr, ranges = r_prime, seqlengths = chrom_lens)
    })
    x_prime <- do.call(c, res)
  } else {
    # here we will send the ranges from multiple
    # chromosomes to a single long chromosome,
    # then perform block bootstrap / permutation,
    # then map back to the original chromosomes.
    stopifnot(all(seg == GenomicRanges::sort(seg)))
    # TODO: code assumes sorted 'x'
    stopifnot(all(x == GenomicRanges::sort(x)))

    new_seg <- map_chroms_to_line_seg(seg)
    new_x <- map_chroms_to_line_seg(x)

    L_c <- sum(chrom_lens)
    seg_length <- seg %>%
      group_by(state) %>%
      summarise(Ls = sum(width)) # derive each states length
    L_s <- seg_length$Ls
    r_prime <- seg_bootstrap_iranges(new_seg, new_x, seg$state, L_c, L_s, L_b, proportion_length, coarse)
    r_prime <- GenomicRanges::sort(r_prime)
    # map the bootstrapped / permuted ranges back to chromosomes:
    x_prime <- map_line_to_chroms_seg(r_prime, x)
  }
  x_prime
}

seg_bootstrap_iranges <- function(seg, x, state, L_c, L_s, L_b, proportion_length = TRUE, coarse = FALSE) {
  if (proportion_length) {
    # number of states
    ns <- length(L_s)

    # the block width for each segmentation state is scaled
    # down to the segmentation state size, e.g. if segmentation state
    # is half of the chromosome, then the block width is half of L_b
    L_b0 <- round(L_b * L_s / L_c)
    obj <- lapply(1:ns, function(m) { # loop over segmentation states
      # the length of the block for this state
      L_bs <- L_b0[m]
      # the segmentation for this state
      seg2 <- seg[state == m]
      # fraction that each range of the segmentation comprises of the whole
      p <- (width(seg2)) / sum(width(seg2))
      # number of blocks within each range of the segmentation
      times <- ceiling(L_c / L_b * p)
      # create random start positions within each segment
      random_start <- unlist(mapply(function(time, x, y) runif(time, x, y), times, start(seg2), end(seg2)))

      # shuffle the blocks
      random_start <- sample(random_start)
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

# these need separate names if they are doing diff things than the ones in 'bootstrap.R'
# (can't have two functions with same name in same R package)
map_chroms_to_line_seg <- function(x) {
  L_c <- seqlengths(x)
  chrom_shift <- c(0, cumsum(as.numeric(unname(L_c)))[-length(L_c)])
  IRanges::shift(ranges(x), rep(chrom_shift, seqnames(x)@lengths))
}

map_line_to_chroms_seg <- function(r_prime, x) {
  L_c <- seqlengths(x)[seqnames(x)@values]
  chrom_shift <- c(0, cumsum(unname(L_c))[-length(L_c)])
  chrom_blocks <- successiveIRanges(width = L_c)
  # r_prime2<-join_overlap_intersect(r_prime,chrom_blocks) # more accurate
  # idx <-findOverlaps(r_prime2,chrom_blocks,select = "first")
  idx <- findOverlaps(r_prime, chrom_blocks, select = "first") # fastest
  r_on_chroms <- IRanges::shift(r_prime, -chrom_shift[idx])
  chroms <- as.character(seqnames(x)@values[idx])
  keep <- start(r_on_chroms) >= 1 & end(r_on_chroms) <= L_c[chroms]
  GRanges(chroms[keep], r_on_chroms[keep], seqlengths = L_c)
  # GRanges(chroms, r_on_chroms, seqlengths = L_c)
}
