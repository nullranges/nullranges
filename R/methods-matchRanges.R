## Internal functions for matchRanges ------------------------------------------------------

#' Validate covariate formula
#' @param f A formula object.
#' @return A validated character vector of covariates
#' @noRd
parseFormula <- function(f) {

  ## Check for valid formula
  if (length(f) != 0) {

    lhs <- rlang::f_lhs(f)
    rhs <- rlang::f_rhs(f)

    if (!is.null(lhs)) {
      msg <- paste0('covar only accepts rhs of formula.')
      stop(msg)
    }

  } else {
    msg <- paste0('covar cannot be an empty formula.')
    stop(msg)
  }

  ## Return covariates
  return(labels(terms(f)))

}

#' Define sample.vec to handle vectors of varying length
#' @noRd
sample.vec <- function(x, ...) x[sample(length(x), ...)]

#' Nearest neighbor matching with replacement (data.table)
#' @param fps A numeric vector of `focal` propensity scores.
#' @param pps A numeric vector of `pool` propensity scores.
#' @inheritParams matchRanges
#' @return A matched `data.table`.
#' @noRd
nnMatch <- function(fps, pps, replace) {

  ## Suppress R CMD CHECK NOTE
  ppsIndex <- NULL
  fpsIndex <- NULL
  N <- NULL

  ## Create data table with original pps index and setkey to sort
  dt <- data.table(pps, val = pps, ppsIndex = seq_along(pps))
  setkey(dt, pps)

  ## Map ids to unique propensity scores
  fpsMap <- data.table(fps = fps, fpsIndex = seq_along(fps))

  ## Order fpsMap by fps (for to adding index later)
  fpsMap <- fpsMap[order(fps)]

  ## Create data table with each unique fps and the number of occurrences
  uniq.fps <- fpsMap[, .(fps, .N), by = fps]

  ## Find all nearest matches for each unique fps
  amdt <- dt[.(uniq.fps), roll = 'nearest', mult = 'all']

  ## Randomly subsample with replacement among duplicates;
  ## Check for categorical variables during sampling
  mdt <-
    amdt[, {
      if (length(N) > 1) {
        stopifnot(length(unique(N)) == 1)
        .(ppsIndex = sample.vec(ppsIndex, N[1], replace = replace))
      } else {
        .(ppsIndex = sample.vec(ppsIndex, N, replace = replace))
      }
    }, by = fps]

  ## Add fpsIndex to mdt
  stopifnot(mdt$fps == fpsMap$fps)
  mdt$fpsIndex <- fpsMap$fpsIndex

  ## Reorder by fpsIndex (to match fps)
  mdt <- mdt[order(fpsIndex)]

  ## Return matched data table
  return(mdt)

}

#' Perform rejection sampling
#' @inheritParams nnMatch
#' @return A binary numeric vector equal to the length
#'   of pps representing whether or not to accept each
#'   pps based on fps distribution.
#' @noRd
rejectSample <- function(fps, pps) {

  ## Ensure fps <= pps
  if (length(fps) > length(pps)) {
    stop("focal must be <= pool for method = 'rejection'.")
  }

  ## Kernal density estimates for focal and pool
  df <- kde(fps)
  dp <- kde(pps)

  ## Set scale by finding the highest point of density ratios (focal/pool)
  ## This ensures that pool covers focal at all points
  grid <- seq(from=quantile(pps, .001), to=quantile(pps, .999), length=1000)
  fgrid <- predict(df, x=grid)
  pgrid <- predict(dp, x=grid)
  if (any(fgrid < 0 | pgrid < 0)) {
    stop("kernel density estimates by ks::kde are negative, cannot perform rejection sampling")
  }
  scale <- max(fgrid/pgrid)
  if (scale > 1e3) {
    stop("scaling factor for density of the PS for pool is > 1e3, could lead to instability")
  }

  ## Calculate the probability of accepting each pool
  thresh <- function(x) ifelse(x > 1e-3, x, 0)
  accept_prob <- pmin(thresh(predict(df, x=pps))/(scale * predict(dp, x=pps)), 1)

  ## Randomly select pps according to accept probability
  accept <- rbinom(length(pps), size=1, prob=accept_prob)

  ## Return random selection equal to the length of pps
  return(accept)

}

#' Rejection sampling and matching for discrete and
#' continuous distributions.
#' @inherit nnMatch params return
#' @noRd
rsMatch <- function(fps, pps, replace) {

  ## Detect descrete distribution
  accept <- tryCatch(expr = {

    ## Rejection sampling
    rejectSample(fps, pps)

  }, error = function(e) {

    ## Introduce random noise to convert discrete into continuous
    stdev <- seq(1e-06, 0.1, length.out = 100)

    ## Iterate over possible noise levels
    for(i in seq_along(stdev)) {

      ## Add noise to fps and pps
      nfps <- fps + (rnorm(length(fps), mean = 0, sd = stdev[i]))
      npps <- pps + (rnorm(length(pps), mean = 0, sd = stdev[i]))

      ## Bound by 0 and 1
      nfps <- pmax(0, pmin(nfps, 1))
      npps <- pmax(0, pmin(npps, 1))

      ## Try rejection sampling
      # Mike: I'm taking off the suppressWarnings,
      #   we can instead catch these upstream?
      accept <- rejectSample(nfps, npps)

      ## Keep trying until there are enough options
      if (any(is.na(accept))) next
      if (sum(accept) >= length(fps)) break

    }

    return(accept)

  })

  ## Ensure there are adequate options
  stopifnot(sum(accept) >= length(fps))

  ## Select matching indices irrespective of data order
  ppsIndex <- sample(which(accept == 1), length(fps), replace = replace)

  ## Assemble matched data table
  mdt <- data.table(fps, ppsIndex, fpsIndex = seq_along(fps))

  ## Return matched data table
  return(mdt)

}

#' Stratifies propensity scores by bins
#' @param fm A `data.table` with `focal` propensity scores
#'   and their corresponding indices.
#' @param fm A `data.table` with `pool` propensity scores
#'   and their corresponding indices.
#' @param n The number of bins used for stratifying
#'   propensity scores.
#' @return A `data.table` of strata with fps, pps, and
#'   their respective indices for each bin.
#' @noRd
stratify <- function(fm, pm, n) {

  fpsIndex <- NULL
  ppsIndex <- NULL
  bin <- NULL

  ## Define breaks using fps and pps
  mn <- min(c(fm$fps, pm$pps))
  mx <- max(c(fm$fps, pm$pps))
  br <- seq(from=mn, to=mx, by=(mx-mn)/n)

  ## Assign fps and pps to bin
  fm$bin <- .bincode(fm$fps, br, include.lowest = TRUE)
  pm$bin <- .bincode(pm$pps, br, include.lowest = TRUE)

  ## Assign indices to bins
  fpsBins <- fm[, .(fpsN = .N, fpsIndices = list(fpsIndex)), by = bin]
  ppsBins <- pm[, .(ppsN = .N, ppsIndices = list(ppsIndex)), by = bin]

  ## Define strata by joining fps and pps on bins
  strata <- fpsBins[ppsBins, on = 'bin']

  return(strata)

}

#' Perform iterative stratified sampling
#' @inherit nnMatch params return
#' @noRd
ssMatch <- function(fps, pps, replace) {
  ## Suppress R CMD CHECK Note
  fpsN <- ppsN <- bin <- NULL 
  fpsIndices <- ppsIndices <- fpsIndex <- ppsIndex <- NULL

  ## Initialize results, fpsOptions and ppsOptions
  results <- data.table(bin=integer(), fpsIndex=integer(), ppsIndex=integer())
  fpsOptions <- data.table(fps, val = fps, fpsIndex = seq_along(fps))
  ppsOptions <- data.table(pps, val = pps, ppsIndex = seq_along(pps))

  ## Set flags and iteration count
  skip <- FALSE
  i <- 1

  ## Start progress bar
  pb <- progress::progress_bar$new(
    format = "  :step [:bar] :percent elapsed: :elapsedfull",
    clear = FALSE, total = length(fps) + 1)
  pb$tick(0)

  while (nrow(results) != length(fps)) {
    ## Update n
    if (skip)  n <- floor(n/2)
    if (!skip) n <- length(unique(c(fpsOptions$fps, ppsOptions$pps)))

    ## Stratify ps by bins and match focal and pool
    strata <- stratify(fpsOptions, ppsOptions, n)
    while (nrow(strata[!is.na(fpsN) & fpsN <= ppsN]) == 0) {
      ## Enter faster n searching
      skip <- TRUE
      n <- floor(n/2)

      ## Stratify ps by bins and match focal and pool
      strata <- stratify(fpsOptions, ppsOptions, n)
    }
    ## Update progress
    pb$update(tokens=list(step=sprintf('Iteration %s, %s bin(s)', i, n)),
              ratio = nrow(results)/length(fps))
    i <- i + 1

    ## Assign indices that can be sampled
    result <-
      strata[!is.na(fpsN) & fpsN <= ppsN,
             .(fpsIndex = unlist(fpsIndices),
               ppsIndex = sample.vec(unlist(ppsIndices), fpsN, replace = replace)),
             by = bin]

    ## Append to results
    results <- rbind(results, result)

    ## Remove assigned indices from options
    fpsOptions <- fpsOptions[!fpsIndex %in% result$fpsIndex]
    ppsOptions <- ppsOptions[!ppsIndex %in% result$ppsIndex]
  }

  ## Close progress bar
  pb$update(tokens=list(step=sprintf('Iteration %s, %s bin(s), done!', i, n)),
            ratio = nrow(results)/length(fps))
  if(pb$finished) pb$terminate()

  ## Reorder by fpsIndex
  results <- results[order(results$fpsIndex)]

  ## Assemble matched data table
  mdt <- data.table(fps = fps[results$fpsIndex],
                    ppsIndex = results$ppsIndex,
                    fpsIndex = results$fpsIndex)

  ## Return matched data table
  return(mdt)
}

#' Calculate propensity scores, implement matching
#' methods, and construct Matched object
#' @param covarData A `data.table` with covariate data.
#'   Contains id column for designating whether covariates
#'   belong in focal or pool (1 or 0, respectively), and
#'   columns for the values of each covariate.
#' @inheritParams matchRanges
#' @return A `Matched` object.
#' @noRd
propensityMatch <- function(covarData, covars, method, replace) {

  ## Assemble covariate formula
  f <- as.formula(paste("id ~", paste(covars, collapse = "+")))

  ## Run glm model

  # April 2023 - speedglm is removed from CRAN so we cannot use it here...
  #model <- speedglm(formula = f, data = covarData,
  #                  family = binomial("logit"), fitted = TRUE, model = TRUE)
  model <- glm(formula = f, data = covarData,
               family = binomial("logit"))
  
  ## Get propensity scores of focal and pool sets as vectors
  id <- NULL
  ps <- NULL

  psData <- data.table(ps = predict(model, type = "response"), id = covarData$id)
  fps <- psData[id == 1, ps]
  pps <- psData[id == 0, ps]

  ## Add propensity scores to covarData
  covarData$ps <- psData$ps

  ## Implement matching method w or w/o replacement
  if (method == 'nearest' & isTRUE(replace))
    mdt <- nnMatch(fps, pps, replace = TRUE)

  if (method == 'nearest' & isFALSE(replace))
    stop("nearest neighbor matching without replacement not available.")

  if (method == 'rejection' & isTRUE(replace))
    mdt <- rsMatch(fps, pps, replace = TRUE)

  if (method == 'rejection' & isFALSE(replace))
    mdt <- rsMatch(fps, pps, replace = FALSE)

  if (method == 'stratified' & isTRUE(replace))
    mdt <- ssMatch(fps, pps, replace = TRUE)

  if (method == 'stratified' & isFALSE(replace))
    mdt <- ssMatch(fps, pps, replace = FALSE)


  ## Assemble information by set
  matchedData <- rbind(
    covarData[id == 1, c(.SD, set = 'focal')],
    covarData[id == 0][mdt$ppsIndex, c(.SD, set = 'matched')],
    covarData[id == 0, c(.SD, set = 'pool')],
    covarData[id == 0][!mdt$ppsIndex, c(.SD, set = 'unmatched')]
  )

  ## Matched indicies
  matchedIndex <- mdt$ppsIndex

  ## Combine information in Matched object
  obj <- Matched(matchedData = matchedData,
                 matchedIndex = matchedIndex,
                 covar = covars,
                 method = method,
                 replace = replace)

  return(obj)
}

#' Core matching for DataFrames/GRanges/GInteractions
#' @param focal A `DataFrame` object containing
#'   the focal data to match.
#' @param pool A `DataFrame` object containing
#'   the pool from which to select matches.
#' @inheritParams matchRanges
#' @return A `Matched` object.
#' @noRd
matchRanges_Core <- function(focal, pool, covar, method, replace) {

  ## Extract covariates from formula as character vector
  covars <- parseFormula(covar)

  ## Check that all covariates are in both focal and pool
  if (!(all(covars %in% colnames(focal)) &
        all(covars %in% colnames(pool)))) {
    stop("All variables in covar must be columns in both focal and pool.")
  }

  ## Check method and replace arguments
  method <- match.arg(method, choices = c('nearest', 'rejection', 'stratified'))

  if (isFALSE(replace) & nrow(focal) >= nrow(pool))
    stop("focal must be <= pool when replace = FALSE.")

  if (method == 'nearest' & isFALSE(replace))
    stop("nearest neighbor matching without replacement not available.")

  if (method == 'stratified' & nrow(focal) >= nrow(pool))
    stop("focal must be <= pool for stratified sampling.")

  ## Create data table with covariate data
  covarData <- as.data.table(cbind(id = factor(c(rep(1, nrow(focal)),
                                                 rep(0, nrow(pool)))),
                                   rbind(focal[covars], pool[covars])))

  ## Calculate propensity scores and match data
  md <- propensityMatch(covarData, covars, method, replace)

  return(md)
}


## Matched subclass methods for matchRanges ----------------------------------------------

#' MatchedDataFrame method for matchRanges
#' @inheritParams matchRanges_Core
#' @return A `MatchedDataFrame` object.
#' @noRd
matchRanges_MatchedDataFrame <- function(focal, pool, covar, method, replace) {

  ## Convert focal and pool to DataFrames
  f <- DataFrame(focal)
  p <- DataFrame(pool)

  ## Execute shared GRanges/GInteractions matching code
  md <- matchRanges_Core(f, p, covar, method, replace)

  ## Combine matched data into MatchedDataFrame class
  obj <- MatchedDataFrame(focal = f,
                          pool  = p,
                          matchedData = matchedData(md),
                          matchedIndex = indices(md),
                          covar = covariates(md),
                          method = method(md),
                          replace = withReplacement(md),
                          delegate = pool[indices(md),])

  ## Return MatchedDataFrame object
  return(obj)
}

#' Generate a covariate-matched control set of ranges
#' 
#' `matchRanges()` uses a propensity score-based method to
#' generate a covariate-matched control set of DataFrame,
#' GRanges, or GInteractions objects.
#' 
#' Available inputs for `focal` and `pool` include `data.frame`,
#' `data.table`, `DataFrame`, `GRanges`, or `GInteractions`.
#' `data.frame` and `data.table` inputs are coerced to `DataFrame`
#' objects and returned as `MatchedDataFrame` while `GRanges` and
#' `GInteractions` objects are returned as `MatchedGRanges` or
#' `MatchedGInteractions`, respectively.
#' 
#' @section Methodology:
#' `matchRanges` uses 
#' [propensity scores](https://en.wikipedia.org/wiki/Propensity_score_matching)
#' to perform subset selection on the `pool` set such that the resulting `matched`
#' set contains similar distributions of covariates to that of the `focal` set.
#' A propensity score is the conditional probability of assigning an element
#' (in our case, a genomic range) to a particular outcome (`Y`) given a set of
#' covariates. Propensity scores are estimated using a logistic regression model
#' where the outcome `Y=1` for `focal` and `Y=0` for `pool`, over the provided 
#' covariates `covar`.
#' 
#' @section Matching methods:
#' \itemize{
#'   \item{`method = 'nearest'`: }{Nearest neighbor matching 
#'   with replacement. Finds the nearest neighbor by using a 
#'   rolling join with `data.table`. Matching without replacement
#'   is not currently supported.}
#'   \item{`method = 'rejection'`: }{(Default) Rejection sampling
#'   with or without replacement. Uses a probability-based approach
#'   to select options in the `pool` that match the `focal` distribition.}
#'   \item{`method = 'stratified'`: }{Iterative stratified sampling
#'   with or without replacement. Bins `focal` and `pool` propensity
#'   scores by value and selects matches within bins until all `focal`
#'   items have a corresponding match in `pool`.}
#' }
#'
#' @param focal A DataFrame, GRanges, or GInteractions object containing
#'   the focal data to match.
#' @param pool A DataFrame, GRanges, or GInteractions object containing
#'   the pool from which to select matches.
#' @param covar A rhs formula with covariates on which to match.
#' @param method A character describing which matching method to use.
#'   supported options are either 'nearest', 'rejection', or 'stratified'.
#' @param replace TRUE/FALSE describing whether to select matches with or without
#'   replacement.
#' @param ... Additional arguments.
#'
#' @return A covariate-matched control set of data.
#'
#' @references
#'
#' matchRanges manuscript:
#'
#' Eric S. Davis, Wancen Mu, Stuart Lee, Mikhail G. Dozmorov,
#' Michael I. Love, Douglas H. Phanstiel. 2022.
#' "matchRanges: Generating null hypothesis genomic ranges
#' via covariate-matched sampling."
#' bioRxiv. doi: 10.1101/2022.08.05.502985
#' 
#' @examples 
#' ## Match with DataFrame
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3)
#' 
#' ## Match with GRanges
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = "GRanges")
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3)
#' 
#' ## Match with GInteractions
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = "GInteractions")
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3)
#' 
#' ## Nearest neighbor matching with replacement
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3,
#'             method = 'nearest',
#'             replace = TRUE)
#' 
#' ## Rejection sampling without replacement
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3,
#'             method = 'rejection',
#'             replace = FALSE)
#' 
#' ## Stratified sampling without replacement
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3,
#'             method = 'stratified',
#'             replace = FALSE)
#'
#' @rdname matchRanges
#'
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @importFrom rlang f_lhs f_rhs
#' @importFrom ks kde
#' @import S4Vectors
#' 
#' @export
setMethod("matchRanges",
          signature = signature(focal   = "DF_OR_df_OR_dt",
                                pool    = "DF_OR_df_OR_dt",
                                covar   = "formula",
                                method  = "character_OR_missing",
                                replace = "logical_OR_missing"),
          definition = matchRanges_MatchedDataFrame)

#' MatchedGRanges method for matchRanges
#' @param focal A `GRanges` object containing
#'   the focal data to match.
#' @param pool A `GRanges` object containing
#'   the pool from which to select matches.
#' @inheritParams matchRanges
#' @return A `MatchedGRanges` object
#' @noRd
matchRanges_MatchedGRanges <- function(focal, pool, covar, method, replace) {

  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)

  ## Execute matching code to get a Matched object
  md <- matchRanges_Core(f, p, covar, method, replace)

  ## Combine matched data into MatchedGRanges class
  obj <- MatchedGRanges(focal = focal,
                        pool = pool,
                        matchedData = matchedData(md),
                        matchedIndex = indices(md),
                        covar = covariates(md),
                        method = method(md),
                        replace = withReplacement(md),
                        delegate = pool[indices(md)])

  ## Return MatchedGRanges object
  return(obj)
}

#' @rdname matchRanges
#' @export
setMethod("matchRanges",
          signature = signature(focal   = "GRanges",
                                pool    = "GRanges",
                                covar   = "formula",
                                method  = 'character_OR_missing',
                                replace = 'logical_OR_missing'),
          definition = matchRanges_MatchedGRanges)

#' MatchedGInteractions method for matchRanges
#' @param focal A `GInteractions` object containing
#'   the focal data to match.
#' @param pool A `GInteractions` object containing
#'   the pool from which to select matches.
#' @inheritParams matchRanges
#' @return A `MatchedGInteractions` object
#' @noRd
matchRanges_MatchedGInteractions <- function(focal, pool, covar, method, replace) {

  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)

  ## Execute matching code to get a Matched object
  md <- matchRanges_Core(f, p, covar, method, replace)

  ## Combine matched data into MatchedGInteractions class
  obj <- MatchedGInteractions(focal = focal,
                              pool = pool,
                              matchedData = matchedData(md),
                              matchedIndex = indices(md),
                              covar = covariates(md),
                              method = method(md),
                              replace = withReplacement(md),
                              delegate = pool[indices(md)])

  ## Return MatchedGInteractions object
  return(obj)
}

#' @rdname matchRanges
#' @export
setMethod("matchRanges",
          signature = signature(focal   = "GInteractions",
                                pool    = "GInteractions",
                                covar   = "formula",
                                method  = 'character_OR_missing',
                                replace = 'logical_OR_missing'),
          definition = matchRanges_MatchedGInteractions)

## Define accessor functions for matchGRanges class --------------------------------------

#' Get focal set from a Matched object
#' 
#' @param x A `MatchedDataFrame`, `MatchedGRanges`,
#'   or `MatchedGInteractions` object.
#' @param ... Additional options.
#' 
#' @return An object of the same class as `x` representing
#'   the focal set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' focal(x)
#' 
#' @rdname focal
#' @export
setMethod("focal", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@focal
})

#' Get pool set from a Matched object
#' 
#' @inheritParams focal
#' 
#' @return An object of the same class as `x` representing
#'   the pool set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' pool(x)
#' 
#' @rdname pool
#' @export
setMethod("pool", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool
})

#' Get matched set from a Matched object
#' 
#' @inheritParams focal
#' 
#' @return An object of the same class as `x` representing
#'   the matched set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' matched(x)
#' 
#' @rdname matched
#' @export
setMethod("matched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@delegate
})

#' Get unmatched set from a Matched object
#' 
#' @inheritParams focal
#' 
#' @return An object of the same class as `x` representing
#'   the unmatched set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' unmatched(x)
#' 
#' @rdname unmatched
#' @export
setMethod("unmatched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool[indices(x, set = "unmatched"),]
})
