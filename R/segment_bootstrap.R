#' Segmentated bootstrap GRanges
#'
#' @param seg segmentation GRanges a column("state") indicate segmentation states
#' @param x the input GRanges
#' @param L_c the length of the chromosome
#' @param L_s a vector of the length of each segmentation region
#' @param L_b the length of the block
#'
#' @export
seg_bootstrap_granges <- function(seg, x, within_chrom = TRUE, L_b) {
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
      x0 <- x[seqnames(x) == chr]
      r_prime <- seg_bootstrap_iranges(seg0, x0, L_c, L_s, L_b)
      GRanges(seqnames = chr, ranges = r_prime, seqlengths = chrom_lens)
    })
    x_prime <- do.call(c, res)
  } else {
    # here we will send the ranges from multiple
    # chromosomes to a single long chromosome,
    # then perform block bootstrap / permutation,
    # then map back to the original chromosomes.

    # TODO: code assumes sorted 'x'
    stopifnot(all(x == GenomicRanges::sort(x)))

    L_c <- sum(chrom_lens)
    seg_length <- seg %>%
      group_by(state) %>%
      summarise(Ls = sum(width)) # derive each states length
    L_s <- seg_length$Ls
    r_prime <- seg_bootstrap_iranges(seg, x, L_c, L_s, L_b)
    r_prime <- GenomicRanges::sort(r_prime)
    # map the bootstrapped / permuted ranges back to chromosomes:
    x_prime <- map_line_to_chroms_seg(r_prime, x)
  }
  x_prime
}

seg_bootstrap_iranges <- function(seg, x, L_c, L_s, L_b) {
  ns <- length(L_s)
  L_b0 <- round(L_s * L_b / L_c) # block width for each segmentation state
  n <- ceiling(L_c / L_b)
  seg <- sort(seg)
  r_seg <- map_chroms_to_line_seg(seg)
  r_x <- map_chroms_to_line_seg(x)

  obj <- lapply(1:ns, function(m) { # loop over segmentation state
    L_bs <- L_b0[m]
    seg2 <- r_seg[seg$state == m]

    p <- seg2@width / sum(seg2@width)
    times <- ceiling(n * p) # number of block within every piece segmentation

    start <- mapply(function(time, x, y) runif(time, x, y), times, seg2@start, end(seg2))
    if (is.list(start)) {
      start <- do.call(c, start)
    }
    start <- sample(start)
    random_blocks <- IRanges(start = start, width = L_bs)

    start_order <- mapply(function(x, y) seq(x, y, L_bs), seg2@start, end(seg2))
    if (is.list(start_order)) {
      start_order <- do.call(c, start_order)
    }
    rearranged_blocks <- IRanges(start = start_order, width = L_bs)

    block_shift <- start(rearranged_blocks) - start(random_blocks)

    fo <- findOverlaps(random_blocks, r_x)
    x_prime <- IRanges::shift(r_x[subjectHits(fo)], block_shift[queryHits(fo)])
    # x_prime<-join_overlap_intersect(x_prime,r_seg) # more accurate
    x_prime <- x_prime[start(x_prime) >= 1 & end(x_prime) <= L_c] # fastest
    x_prime
  })
  obj_prime <- do.call(c, obj)
  return(obj_prime)
}

# these need separate names if they are doing diff things than the ones in 'bootstrap.R'
# (can't have two functions with same name in same R package)
map_chroms_to_line_seg <- function(x) {
  L_c <- seqlengths(x)[seqnames(x)@values]
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
