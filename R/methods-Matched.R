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

#' Internal function for parsing common plotting args
#' 
#' Defines colors & linetypes, parses the set argument
#' and extracts matched data.
#' 
#' @inheritParams plotCovariate
#' @return A list of arguments:
#' * [`md`] - matched data
#' * [`cols`] - named color vector
#' * [`sets`] - parsed sets
#' * [`lty`] = named linetype vector
#' @noRd
parse_plot_args <- function(sets, x){
  
  ## Define colors & linetypes
  cols <- c(focal = "#1F78B4", matched = "#A6CEE3",
            pool = "#33A02C", unmatched = "#B2DF8A")
  lty <- c(focal = 1, matched = 2, pool = 1, unmatched = 2)
  
  ## Parse arguments
  sets <- match.arg(sets, choices = names(cols), several.ok = TRUE)
  
  ## Extract matched data & subset by sets
  md <- matchedData(x)
  md <- md[md[["set"]] %in% sets]
  cols <- cols[names(cols) %in% sets]
  
  ## Return args
  return(list(md = md, cols = cols, sets = sets, lty = lty))
}

#' Internal function for setting type argument
#' @inheritParams plotCovariate
#' @param data Matched data
#' @param x Character or symbol describing which column of matched data
#'   for which to apply function.
#' @return A character with the new value for `type`
#' @noRd
set_matched_type <- function(data, type, thresh, x) {
  
  ## Extract covariate values for checking data type
  vals <- data[[x]]
  
  ## Use class of covariate to set plot type
  if (is.null(type)) {
    if (is.numeric(vals)) {
      type <- ifelse(uniqueN(vals) <= thresh, 'bars', 'lines')
    } else if (is.factor(vals) |
               is.logical(vals) |
               is.character(vals)) {
      type <- 'bars'
    } else {
      stop(sprintf('Data type %s is not supported.'), vals)
    }
  }
  
  type
}

#' Internal function for setting matched plot
#' @inheritParams plotCovariate
#' @param data Matched data
#' @param cols,lty Named character vector of colors/lty for focal,
#'   matched, pool, and unmatched sets.
#' @param x Character or symbol describing which column of matched data
#'   for which to apply function.
#' @importFrom rlang !! enquo
#' @return A `ggplot` set by `type` argument
#' @noRd
set_matched_plot <- function(data, type, cols, lty, thresh, x) {
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
    ans <- ggplot(data, mapping = aes(x = !!x, color = !!set, linetype=!!set)) +
      stat_density(geom = 'line', position = 'identity', na.rm = TRUE) +
      scale_color_manual(values = cols) +
      scale_linetype_manual(values = lty)
  }
  
  if (identical(type, "bars")) {
    ## Suppress R CMD CHECK Note
    N <- pct <- V1 <- NULL
    
    ## Convert categorical-numeric to factor
    vals <- data[[x]]
    if (is.numeric(vals) & uniqueN(vals) <= thresh){
      data[[x]] <- as.factor(signif(vals, 2))
    }
    
    ## Form melted table, calculate percentages, and order (for continuous)
    data <- data[, .N, by = .(eval(set), eval(x))]
    data <- data[, .(eval(x), "pct" = (N/sum(N)*100)), by = set]
    data <- data[order(V1)]
    
    ## Rename columns
    colnames(data) <- c(deparse(set), deparse(x), 'pct')
    
    ## Plot
    ans <- ggplot(data = data, mapping = aes(x = !!set, y = pct, fill = !!x)) +
      geom_col(position = 'stack') +
      labs(y = "Percentage") + 
      scale_fill_hue(l = 70, c = 50)
  }
  
  ## Apply general plot formatting
  ans <- ans +
    theme_minimal()+
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = 'transparent'))
  
  ans
}

#' Internal function to apply log-transformation
#' @inheritParams plotCovariate
#' @param ans Input `ggplot` to be log-transformed
#' @param x Character or symbol describing which column of matched data
#'   for which to apply function.
#' @return A log-transformed `ggplot`.
#' @noRd
apply_log_trans <- function(log, type, ans, x) {
  
  ## Parse log parameter
  if (!is.null(log)) {
    log <- match.arg(log, choices = c('x', 'y'), several.ok = TRUE)
  }
  
  ## Apply x-transformation
  if ('x' %in% log) {
    if (identical(type, 'bars')) {
      stop("Transformation of x-axis not valid when type = 'bars'.")
    }
    ans <- ans +
      scale_x_continuous(
        trans = "log",
        breaks = scales::log_breaks(base = exp(1)),
        oob = scales::oob_squish_infinite) +
      labs(x = sprintf("log(%s)", x))
  }
  
  ## Apply y-transformation
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

#' Internal function for plotting propensity scores
#' @noRd
plot_propensity <- function(x, sets, type, log, thresh = 12) {
  
  ## Suppress R CMD Check Note
  md <- cols <- lty <- NULL

  ## Extract matchedData (md); parse cols, lty & sets
  list2env(parse_plot_args(sets, x), envir = environment())
  
  ## Use class of covariate to set type
  type <- set_matched_type(data = md,
                           type = type,
                           thresh = thresh,
                           x = "ps")
  
  ## Generate plot by type
  ans <- set_matched_plot(data = md,
                          type = type,
                          cols = cols,
                          lty = lty,
                          thresh = thresh,
                          x = "ps")
  
  ## Apply general plot formatting
  ans <- ans +
    labs(title = paste0("~", paste(covariates(x), collapse = "+")))
  
  if (!identical(type, 'bars')) {
    ans <- ans +
      labs(x = "Propensity Score", y = NULL)  
  }
  
  
  ## Apply any log-transformations
  ans <- apply_log_trans(log = log,
                         type = type,
                         ans = ans,
                         x = "ps")
  
  ans

}

#' Internal function for plotting a covariate
#' @noRd
plot_covariate <- function(x, covar, sets, type, log, thresh = 12) {
  
  ## Suppress R CMD Check Note
  md <- cols <- lty <- NULL
  
  ## Extract matchedData (md); parse cols, lty & sets
  list2env(parse_plot_args(sets, x), envir = environment())
  
  ## Parse covariate arguments
  covar <- match.arg(covar, choices = covariates(x), several.ok = FALSE)
  
  ## Use class of covariate to set type
  type <- set_matched_type(data = md,
                           type = type,
                           thresh = thresh,
                           x = covar)
  
  ## Generate plot by type
  ans <- set_matched_plot(data = md,
                          type = type,
                          cols = cols,
                          lty = lty,
                          thresh = thresh,
                          x = !!covar)
  
  ## Apply any log-transformations
  ans <- apply_log_trans(log = log,
                         type = type,
                         ans = ans,
                         x = covar)
  
  ans
  
}

#' Propensity score plotting for Matched objects
#'
#' This function plots the distribution of propensity scores
#' from each matched set of a Matched object.
#' 
#' `plotPropensity` uses the `thresh` argument
#' to determine whether to plot propensity scores as 
#' continuous (line plot) or catetgorical (bar plot).
#' These settings can also be overwritten manually.
#' 
#' @inheritParams plotCovariate
#' 
#' @return Returns a plot of propensity score distributions
#' among matched sets.
#' 
#' @examples
#' ## Matched example dataset
#' mdf <- makeExampleMatchedDataSet(matched = TRUE)
#' 
#' ## Visualize propensity scores
#' plotPropensity(mdf)
#' plotPropensity(mdf,
#'               sets = c('focal', 'matched', 'pool'))
#' plotPropensity(mdf,
#'               sets = c('focal', 'matched', 'pool'),
#'               type = 'ridges')
#' plotPropensity(mdf,
#'               sets = c('focal', 'matched', 'pool'),
#'               type = 'jitter')
#'
#' @rdname plotPropensity
#' @seealso [plotCovariate()] to plot covariate distributions.
#' @import ggplot2 ggridges
#' @export
setMethod("plotPropensity",
          signature = signature(x="Matched",
                                sets = 'character_OR_missing',
                                type = 'character_OR_missing',
                                log = 'character_OR_missing'),
          definition = plot_propensity)

#' Covariate plotting for Matched objects
#' 
#' This function plots the distributions of a covariate
#' from each matched set of a Matched object.
#' 
#' By default, `plotCovariate` will sense the 
#' class of covariate and make a plot best suited to
#' that data type. For example, if the covariate class
#' is categorical in nature then the `type` argument
#' defaults to 'bars'. `type` is set to 'lines' for
#' continuous covariates. These settings can also be overwritten
#' manually.
#' 
#' @param x Matched object
#' @param covar Character naming the covariate to plot.
#'   If multiple are provided, only the first one is used.
#' @param sets Character vector describing which matched set(s)
#'   to include in the plot. Options are 'focal', 'matched',
#'   'pool', or 'unmatched'. Multiple options are accepted.
#' @param type Character naming the plot type. Available
#'   options are one of either 'ridges', 'jitter', 'lines',
#'   or 'bars'. Note that for large datasets, use of 'jitter'
#'   is discouraged because the large density of points can
#'   stall the R-graphics device.
#' @param log Character vector describing which axis or
#'   axes to apply log-transformation. Available options are
#'   'x' and/or 'y'.
#' @param thresh Integer describing the number of
#'   unique values required to classify a numeric
#'   variable as discrete (and convert it to a factor).
#'   If the number of unique values exceeds `thresh`
#'   then the variable is considered continuous.
#' @param ... Additional arguments. 
#' 
#' @return Returns a plot of a covariate's distribution
#'   among matched sets.
#' 
#' @examples
#' ## Matched example dataset
#' mdf <- makeExampleMatchedDataSet(matched = TRUE)
#' 
#' ## Visualize covariates
#' plotCovariate(mdf)
#' plotCovariate(mdf, covar = 'covar2')
#' plotCovariate(mdf,
#'               covar = 'covar1',
#'               sets = c('focal', 'matched', 'pool'))
#' plotCovariate(mdf,
#'               covar = 'covar1',
#'               sets = c('focal', 'matched', 'pool'),
#'               type = 'ridges')
#' plotCovariate(mdf,
#'               covar = 'covar1',
#'               sets = c('focal', 'matched', 'pool'),
#'               type = 'jitter')
#'
#' @rdname plotCovariate
#' @seealso [plotPropensity()] to plot propensity scores.
#' @family Matched plotting
#' @importFrom scales squish_infinite
#' @export
setMethod("plotCovariate",
          signature = signature(x="Matched",
                                covar = 'character_OR_missing',
                                sets = 'character_OR_missing',
                                type = 'character_OR_missing',
                                log = 'character_OR_missing'),
          definition = plot_covariate)
