## Define accessor functions for Matched class -------------------------------------------

#' Accessor methods for Matched Class
#'
#' Functions that get data from Matched subclasses
#' such as Matched, MatchedDataFrame, MatchedGRanges,
#' and MatchedGInteractions.
#'
#' @param x Matched object
#' @param ... additional arguments
#'
#' @name Matched
#' @rdname matched
NULL

#' @rdname matched
#' @export
setMethod("matchedData", "Matched", function(x, ...) {
  x@matchedData
})

#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname matched
#' @export
setMethod("covariates", "Matched", function(x, ...) {
  x@covar
})

getIndices <- function(x, group = 'matched') {
  ## Get group argument
  group <- match.arg(group, choices=c("focal","matched","pool","unmatched"))

  ## Get the length of each focal and pool
  n.focal <- nrow(x@matchedData[group == 'focal'])
  n.pool  <- nrow(x@matchedData[group == 'pool'])

  if (group == 'focal')
    return(1:n.focal)

  if (group == 'matched')
    return(x@matchedIndex)

  if (group == 'unmatched')
    return(which(!(1:n.pool) %in% x@matchedIndex))

  if (group == 'pool')
    return(1:n.pool)
}

#' @param group a character string describing from which group to extract indices.
#'              can be one of 'focal', 'matched', 'pool', or 'unmatched'.
#' @rdname matched
#' @export
setMethod("indices", "Matched", getIndices)


## Define overview method for Matched class ----------------------------------------------

overviewMatched <- function(x) {

  ## Define aggregation function
  agg <- function(x) {
    if (any(is.factor(x)))
      list(table(x))
    else
      list(mean = mean(x), sd = sd(x))
  }


  ## Apply aggregation to matchedData
  md.agg <- x@matchedData[, as.list(c(N =.N, unlist(lapply(.SD, agg)))),
                          .SDcols = -c('id'), by = group]

  ## Calculate distances between focal and matched
  d <- x@matchedData[group == 'focal', -c('id', 'group')] -
    x@matchedData[group == 'matched', -c('id', 'group')]

  ## Apply aggregation to distances
  d.agg <- d[, as.list(unlist(lapply(.SD, agg)))]

  ## Display overview
  cat(class(x), "object:", '\n', sep = ' ')
  print(md.agg, row.names = FALSE)
  cat('--------\n')
  cat('focal - matched: \n')
  print(d.agg, row.names = FALSE)

}

#' @rdname matched
#' @export
setMethod("overview", signature(x="Matched"), overviewMatched)


## Define plot methods for Matched class -------------------------------------------------

#' @importFrom rlang !! enquo
# internal function for covariate plotting
set_matched_plot <- function(type, data, cols, x) {
  x <- enquo(x)
  switch(type,
         "jitter" =
           ggplot(data, aes = aes(x = !!x, y = group, color = group)) +
            geom_jitter(height = 0.25, width = 0, alpha = 0.7) +
            scale_color_manual(values = cols),
         "ridge" =
           ggplot(data,  aes = aes(x = !!x, y = group, fill = group)) +
            geom_density_ridges(alpha = 0.7, color = NA) +
            scale_fill_manual(values = cols),
         "lines" =
           ggplot(data, aes = aes(x = !!x, color = group)) +
            geom_density(show.legend = FALSE, na.rm = TRUE) +
            stat_density(geom = 'line', position = 'identity', na.rm = TRUE) +
            scale_color_manual(values = cols)
  )

}

## Define function for plotting propensity scores
plot_propensity <- function(x, type = NULL) {

  ## Extract matchedData
  md <- matchedData(x)

  ## Define colors
  cols <- c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A")

  if (is.null(type)) {
    type <- ifelse(nrow(md[group == "pool"]) <= 10000, "jitter", "ridge")
  }

  ans <-set_matched_plot(md, type, cols, x = ps)

  ans +
    scale_y_discrete(limits = rev) +
    labs(x = "Propensity Score", y = NULL) +
    theme_minimal() +
    theme(legend.position = 'none',
          panel.border = element_rect(fill = 'transparent'))


}

## Define function for plotting covariates
plot_covariates <- function(x, covar = 'all', sets = 'all', type = NULL, logTransform = FALSE) {

  ## Extract matchedData
  md <- matchedData(x)

  ## Chose plot type by data size
  if (is.null(type)) {
    type <- ifelse(nrow(md[group == "pool"]) <= 10000, "jitter", "ridge")
  }

  ## Define colors & covariates
  cols <- c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A")
  names(cols) <- c('focal', 'matched', 'pool', 'unmatched')

  ## Parse covariate & set to plot
  covar <- match.arg(covar, choices = c('all', covariates(x)), several.ok = TRUE)
  sets <- match.arg(sets, choices = c('all', names(cols)), several.ok = TRUE)

  if (length(covar) == 1) {
    if (covar == 'all') {
      covar <- covariates(x)
    }
  }

  if (length(sets) == 1) {
    if (sets == 'all') {
      sets <- names(cols)
    }
  }

  ## Melt data for plotting multiple covariates
  mmd <- melt(md, measure.vars = covar)

  ## Subset group by sets
  mmd <- mmd[group %in% sets]

  # TODO this leads to 'no visible binding for global variable'
  # produces NOTE during package check
  ans <- set_matched_plot(data,
                          type,
                          cols = cols[names(cols) %in% sets],
                          x = value)

  ans <- ans +
    scale_y_discrete(limits = rev) +
    labs(y = "") +
    facet_grid(~variable, scales = "free_x")+
    theme_minimal()+
    theme(legend.position = 'none',
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = 'transparent'))

  if (logTransform) {
    ans <- ans +
      scale_x_continuous(
        trans = "log",
        breaks = scales::log_breaks(base = exp(1)),
        oob = scales::oob_squish_infinite) +
      labs(x = "log(value)",  y = "")

  }

  ans

}

#' @title Plotting functions for Matched objects
#'
#' @rdname matched-plotting
#' @importFrom graphics plot
#' @import ggplot2 ggridges
#' @method plot Matched
#' @export
setMethod("plot", signature(x="Matched", y="missing"), plot_propensity)

#' @param x ...
#' @param covar ...
#' @param sets ...
#' @param type ...
#' @param logTransform ...
#' @param ... additional arguments
#'
#' @rdname matched-plotting
#' @importFrom scales squish_infinite
#' @export
setMethod("plotCovariates", signature(x="Matched"), function(x, ...) {
  # not sure why you need NSE here?
  mc <- match.call()
  mc[[1]] <- quote(plot_covariates)
  eval(mc)
})
