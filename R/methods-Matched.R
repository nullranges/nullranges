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

getIndices <- function(x, set = 'matched') {
  ## Get set argument
  set <- match.arg(set, choices=c("focal","matched","pool","unmatched"))

  ## Get the length of each focal and pool
  n.focal <- nrow(x@matchedData[set == 'focal'])
  n.pool  <- nrow(x@matchedData[set == 'pool'])

  if (set == 'focal')
    return(1:n.focal)

  if (set == 'matched')
    return(x@matchedIndex)

  if (set == 'unmatched')
    return(which(!(1:n.pool) %in% x@matchedIndex))

  if (set == 'pool')
    return(1:n.pool)
}

#' @param set a character string describing from which set to extract indices.
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
  set <- NULL

  md <- matchedData(x)
  ## Apply aggregation to matchedData
  md.agg <- md[, as.list(c(N =.N, unlist(lapply(.SD, agg)))),
                          .SDcols = -c('id'), by = set]

  ## Calculate distances between focal and matched
  d <- md[set == 'focal', -c('id', 'set')] -
    md[set == 'matched', -c('id', 'set')]

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
set_matched_plot <- function(data, type, cols, x) {
  x <- rlang::ensym(x)
  set <- rlang::sym("set")
  
  type <- match.arg(type, c("jitter", "ridges", "lines", "bars"))
  
  if (identical(type, "jitter")) {
    ans <- ggplot(data, mapping = aes(x = !!x, y = !!set, color = !!set)) +
      geom_jitter(height = 0.25, width = 0, alpha = 0.7) +
      scale_y_discrete(limits = rev) +
      scale_color_manual(values = cols)
  }
  
  if (identical(type, "ridges")) {
    ans <-
      ggplot(data,  mapping = aes(x = !!x, y = !!set, fill = !!set)) +
      geom_density_ridges(alpha = 0.7, color = NA) +
      scale_y_discrete(limits = rev) +
      scale_fill_manual(values = cols)
  }
  
  if (identical(type, "lines")) {
    ans <- ggplot(data, mapping = aes(x = !!x, color = !!set)) +
      geom_density(show.legend = FALSE, na.rm = TRUE) +
      stat_density(geom = 'line', position = 'identity', na.rm = TRUE) +
      scale_color_manual(values = cols)
  }
  
  if (identical(type, "bars")) {
    ## Form melted table, calculate percentages, and order (for continuous)
    data <- data[, .N, by = .(eval(set), eval(x))]
    data <- data[, .(eval(x), pct = (N/sum(N)*100)), by = set]
    data <- data[order(V1)]
    
    ## Rename columns
    colnames(data) <- c(deparse(set), deparse(x), 'pct')
    
    ## Plot
    ans <- ggplot(data = data, mapping = aes(x = !!set, y = pct, fill = !!x)) +
      geom_col(position = 'stack') +
      labs(y = "Percentage")
  }
  
  ans
}

## Define function for plotting propensity scores
plot_propensity <- function(x, type = NULL) {

  ## Extract matchedData
  md <- matchedData(x)

  ## Define colors
  cols <- c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A")



  if (is.null(type)) {
    type <- ifelse(sum(md[["set"]] == "pool") <= 10000, "jitter", "ridge")
  }

  ans <- set_matched_plot(md, type, cols, x = "ps")

  ans +
    labs(x = "Propensity Score", y = NULL) +
    theme_minimal() +
    theme(legend.position = 'none',
          panel.border = element_rect(fill = 'transparent'))


}

## Define function for plotting covariates
plot_covariate <- function(x, covar = NULL, sets = 'all', type = NULL, log = NULL) {

  ## Define colors & sets
  cols <- c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A")
  names(cols) <- c('focal', 'matched', 'pool', 'unmatched')
  
  ## Parse arguments
  covar <- match.arg(covar, choices = covariates(x), several.ok = FALSE)
  sets <- match.arg(sets, choices = c('all', names(cols)), several.ok = TRUE)
  
  if (length(sets) == 1) {
    if (sets == 'all') {
      sets <- names(cols)
    }
  }
  
  if (!is.null(log)) {
    log <- match.arg(log, choices = c('x', 'y'), several.ok = TRUE)
  }
  
  ## Extract matched data
  md <- matchedData(x)
  
  ## Subset matched data and colors by sets
  md <- md[md[["set"]] %in% sets]
  cols <- cols[names(cols) %in% sets]
  
  ## Extract covariate values for checking data type
  covarValues <- md[[covar]]
  
  ## Use class of covariate values to set plot type
  if (is.null(type)) {
    if (is.numeric(covarValues)) {
      type <- ifelse(uniqueN(covarValues) <= 10, 'bars', 'lines')
    } else if (is.factor(covarValues) |
               is.logical(covarValues) |
               is.character(covarValues)) {
      type <- 'bars'
    } else {
      stop(sprintf('Data type %s is not supported.'), covarValues)
    }
  }
  
  ## Generate plot by type
  ans <- set_matched_plot(data = md,
                          type = type,
                          cols = cols,
                          x = !!covar)
  
  ## Apply general plot formatting
  ans <- ans +
    theme_minimal()+
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = 'transparent'))
  
  ## Apply log transformation(s)
  if ('x' %in% log) {
    if (identical(type, 'bars')) {
      stop("Transformation of x-axis not valid when type = 'bars'.")
    }
    ans <- ans +
      scale_x_continuous(
        trans = "log",
        breaks = scales::log_breaks(base = exp(1)),
        oob = scales::oob_squish_infinite) +
      labs(x = sprintf("log(%s)", covar))
  }
  
  if ('y' %in% log) {
    if (!identical(type, 'lines')) {
      stop("Transformation of y-axis only valid when type = 'lines'.")
    }
    ans <- ans +
      scale_y_continuous(
        trans = "log",
        breaks = scales::log_breaks(base = exp(1)),
        oob = scales::oob_squish_infinite) +
      labs(y = "log(density)",  y = "")
  }
  
  ans

}

#' @title Plotting functions for Matched objects
#'
#' @rdname matched-plotting
#' @import ggplot2 ggridges
#' @export
setMethod("plotPropensity", signature(x="Matched"),
          function(x, type = NULL) plot_propensity(x, type))

#' @param x ...
#' @param covar ...
#' @param sets ...
#' @param type ...
#' @param log ...
#' @param ... additional arguments
#'
#' @rdname matched-plotting
#' @importFrom scales squish_infinite
#' @export
setMethod("plotCovariate", signature(x="Matched"), function(x, covar = NULL, sets = 'all', type = NULL, log = NULL) {
  # not sure why you need NSE here?
  plot_covariate(x, covar, sets, type, log)
})
