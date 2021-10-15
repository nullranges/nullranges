#' Genomic segmentation based on gene density
#'
#' @param seg the segmentation GRanges returned by \code{segmentDensity} function
#' @param exclude GRanges of excluded region
#' @param type the type of plot returned. Choices are segmentation plot
#' included ranges information, barplot showing segmentation states' distribution across chromosome,
#' or a box plot indicating average density within each states.
#' Default is all plots are displayed.
#' The y axis \code{"density"} represent square root of overlap counts within segment length. If a user defined
#' segmentation GRanges is given, the y axis is default to be 1 for all states in segmentation plot.
#'  
#' @param region GRanges of stricted region that want to be plotted.
#'
#' @return A `ggplot` set by `type` argument
#'
#' @importFrom IRanges findOverlaps
#' @importFrom plyranges summarise group_by %>%
#' @importFrom scales breaks_extended label_comma
#' @importFrom grDevices palette.colors
#'
#' @examples
#'
#' example("segmentDensity")
#' plotSegment(seg, exclude, type = "ranges")
#' plotSegment(seg, exclude, type = "barplot")
#' plotSegment(seg, exclude, type = "boxplot")
#' @export
plotSegment <- function(seg, exclude = NULL,
                        type = c("ranges", "barplot", "boxplot"),
                        region = NULL) {
  type <- match.arg(type, c("ranges", "barplot", "boxplot"))

  if ("counts" %in% colnames(mcols(seg))) {
    counts <- mcols(seg)$counts
  } else {
    counts <- 1
  }
  
  if (!is.null(exclude)) {
    mcols(exclude)$state <- "excluded"
  }

  if (!is.null(region)) {
    seg_fo <- findOverlaps(seg, region)
    seg <- join_overlap_intersect(seg[queryHits(seg_fo)], region)
    counts <- counts[queryHits(seg_fo)]
    if (!is.null(exclude)) {
      exclude <- join_overlap_intersect(exclude, region)
    }
  }

  if (!is.null(exclude)) {
    full_query <- c(seg, exclude)
    len_exclude <- length(exclude)
  } else {
    full_query <- seg
    len_exclude <- 0
  }
  
  q <- quantile(sqrt(counts), .975)
  seq2 <- pmin(sqrt(counts), q)

  dat <- data.frame(
    chr = seqnames(full_query),
    counts = c(seq2, rep(0, len_exclude)),
    state = factor(full_query$state),
    start = IRanges::start(full_query),
    end = IRanges::end(full_query),
    width = IRanges::width(full_query)
  )

  cols <- unname(palette.colors()[-c(1,5)])

  if (identical(type, "ranges")) {
    ans <- ggplot2::ggplot() +
      labs(col = "States", y = "density") +
      geom_segment(aes_string(
        x = "start", y = "counts",
        xend = "end", yend = "counts", col = "state"
      ),
      size = 2.5, data = dat, show.legend = TRUE
      ) +
      scale_color_manual(values = cols) +
      facet_wrap(~chr, ncol = 6) +
      scale_y_continuous(breaks = scales::breaks_extended(4))
    if (!is.null(region)) {
      ans <- ans + scale_x_continuous(
        name = "Location(b)",
        labels = label_comma()
      )
    } else {
      ans <- ans + scale_x_continuous(
        name = "Location(Mb)",
        breaks = c(0.5e8, 1.5e8, 2.5e8),
        labels = c("50", "150", "250")
      )
    }
  }

  if (identical(type, "barplot")) {
    dat_bar <- dat %>%
      plyranges::group_by(.data$chr, .data$state) %>%
      plyranges::summarise(distribution = sum(.data$width))
    ans <- ggplot(dat_bar, aes(fill = .data$state,
                               y = .data$distribution, x = .data$chr)) +
      geom_bar(position = "fill", stat = "identity") +
      theme_classic() +
      scale_fill_manual(values = cols) +
      scale_x_discrete(guide = guide_axis(angle = 45))
  }

  if (identical(type, "boxplot")) {
    states <- data.frame(state = factor(seg$state), count = seq2)
    ans <- ggplot2::ggplot(aes_string(x = "state", y = "count", fill = "state"),
                           data = states) +
      geom_boxplot() +
      scale_fill_manual(values = cols) +
      labs(x = "States", y = "density")
  }
  ## Apply general plot formatting
  ans <- ans +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = "transparent")
    )

  ans
}

#' Combine nearby regions with same state
#'
#' @param x the input GRanges
#' @param col the name of the column for the segment states
#'
#' @return a GRanges with metadata column \code{state}
#' giving the segmentation state
#'
#' @examples
#'
#' n <- 10000
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges(round(
#'   c(runif(n/4,1,991), runif(n/4,1001,3991),
#'     runif(n/4,4001,4991), runif(n/4,7001,9991))),
#'   width=10), seqlengths=c(chr1=10000))
#' gr$name <- rep(1:4,each=10)
#' gr <- sort(gr)
#' seg <- reduceSegment(gr, col="name")
#'
#' @export
reduceSegment <- function(x, col="state") {

  n <- max(mcols(x)[,col])
  
  seg <- do.call(c, lapply(seq_len(n), function(s) {
    x <- reduce(x[mcols(x)[,col] == s])
    mcols(x)$state <- s
    x
  }))

  seg <- sortSeqlevels(seg)
  seg <- GenomicRanges::sort(seg)

  return(seg)
}
