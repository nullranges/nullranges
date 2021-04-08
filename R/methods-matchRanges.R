## Helper functions for matchRanges ------------------------------------------------------

## Define helper function to parse formula
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

## Define sample.vec to handle vectors of varying length
sample.vec <- function(x, ...) x[sample(length(x), ...)]

## Helper function that calculates propensity scores
## and implements nearest neighbor matching
propensityMatch <- function(covarData, covars) {
  
  ## Assemble covariate formula
  f <- as.formula(paste("id ~", paste(covars, collapse = "+")))
  
  ## Run glm model
  model <- glm(formula = f, data = covarData, family = binomial("logit"))
  
  ## Get propensity scores of focal and pool groups as vectors
  psData <- data.table(ps = predict(model, type = "response"), id = model$model$id)
  fps <- psData[id == 1, ps]
  pps <- psData[id == 0, ps]
  
  ## Add propensity scores to covarData
  covarData$ps <- psData$ps
  
  ## Create data table with original pps index and setkey to sort
  dt <- data.table(pps, val = pps, ppsIndex = 1:length(pps))
  setkey(dt, pps)
  
  ## Map ids to unique propensity scores
  fpsMap <- data.table(fps = fps, fpsIndex = seq_along(fps))

  ## Order fpsMap by fps (for to adding index later)
  fpsMap <- fpsMap[order(fps)]

  ## Create data table with each unique fps and the number of occurrences
  uniq.fps <- fpsMap[, .(fps, .N), by = fps]

  ## Find all nearest matches for each unique fps
  amdt <- dt[.(uniq.fps), roll = 'nearest', mult = 'all']

  ## Randomly subsample with replacement among duplicates
  set.seed(123)
  mdt <- amdt[, .(ppsIndex = sample.vec(ppsIndex, N, replace = T)), by = fps]

  ## Add fpsIndex to mdt
  stopifnot(mdt$fps == fpsMap$fps)
  mdt$fpsIndex <- fpsMap$fpsIndex

  ## Reorder by fpsIndex (to match fps)
  mdt <- mdt[order(fpsIndex)]

  ## Assemble information by group
  matchedData <- rbind(
    covarData[id == 1, c(.SD, group = 'focal')],
    covarData[id == 0][mdt$ppsIndex, c(.SD, group = 'matched')],
    covarData[id == 0, c(.SD, group = 'pool')],
    covarData[id == 0][!mdt$ppsIndex, c(.SD, group = 'unmatched')]
  )
  
  ## Matched indicies
  matchedIndex <- mdt$ppsIndex
  
  ## Combine information in Matched object
  obj <- Matched(matchedData = matchedData,
                 matchedIndex = matchedIndex,
                 covar = covars)
  
  return(obj)
}

## Helper function - matching core for DataFrames/GRanges/GInteractions
matchRanges.Core <- function(focal, pool, covar) {
  
  ## Extract covariates from formula as character vector
  covars <- parseFormula(covar)
  
  ## Check that all covariates are in both focal and pool
  if (!(all(covars %in% colnames(focal)) &
        all(covars %in% colnames(pool)))) {
    stop("All variables in covar must be columns in both focal and pool.")
  }
  
  ## Create data table with covariate data
  covarData <- as.data.table(cbind(id = factor(c(rep(1, nrow(focal)),
                                                 rep(0, nrow(pool)))),
                                   rbind(focal[covars], pool[covars])))
  
  ## Calculate propensity scores and match data
  md <- propensityMatch(covarData, covars)
  
  return(md)
}


## Matched subclass methods for matchRanges ----------------------------------------------

## MatchedDataFrame method for matchRanges
matchRanges.MatchedDataFrame <- function(focal, pool, covar) {
  
  ## Convert focal and pool to DataFrames
  f <- DataFrame(focal)
  p <- DataFrame(pool)
  
  ## Execute shared GRanges/GInteractions matching code
  md <- matchRanges.Core(f, p, covar)
  
  ## Combine matched data into MatchedDataFrame class
  obj <- MatchedDataFrame(focal = f,
                          pool  = p,
                          matchedData = matchedData(md),
                          matchedIndex = indices(md),
                          covar = covariates(md),
                          delegate = pool[indices(md),])
  
  ## Return MatchedDataFrame object
  return(obj)
}

#' matchRanges
#'
#' @description 
#' A function that generates a covariate-matched control
#' set of the same datatype.
#'
#' @param focal a DataFrame, GRanges, or GInteractions object containing
#'              the focal data to match.
#' @param pool  a DataFrame, GRanges, or GInteractions object containing
#'              the pool from which to select matches.
#' @param covar a rhs formula with covariates on which to match.
#' @param ...   additional arguments
#' 
#' @return a covariate-matched control set of data
#'
#' @rdname matchRanges
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @importFrom rlang f_lhs
#' @importFrom rlang f_rhs
#' @import S4Vectors
#' @export
setMethod("matchRanges",
          signature = signature(focal = "DF_OR_df_OR_dt",
                                pool  = "DF_OR_df_OR_dt",
                                covar = "formula"),
          definition = matchRanges.MatchedDataFrame)


## MatchedGRanges method for matchRanges
matchRanges.MatchedGRanges <- function(focal, pool, covar) {
  
  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)
  
  ## Execute matching code to get a Matched object
  md <- matchRanges.Core(f, p, covar)
  
  ## Combine matched data into MatchedGRanges class
  obj <- MatchedGRanges(focal = focal,
                        pool = pool,
                        matchedData = matchedData(md),
                        matchedIndex = indices(md),
                        covar = covariates(md),
                        delegate = pool[indices(md)])
  
  ## Return MatchedGRanges object
  return(obj)
}

#' @rdname matchRanges
#' @export
setMethod("matchRanges",
          signature = signature(focal = "GRanges",
                                pool  = "GRanges",
                                covar = "formula"),
          definition = matchRanges.MatchedGRanges)


## MatchedGInteractions method for matchRanges
matchRanges.MatchedGInteractions <- function(focal, pool, covar) {
  
  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)
  
  ## Execute matching code to get a Matched object
  md <- matchRanges.Core(f, p, covar)
  
  ## Combine matched data into MatchedGInteractions class
  obj <- MatchedGInteractions(focal = focal,
                              pool = pool,
                              matchedData = matchedData(md),
                              matchedIndex = indices(md),
                              covar = covariates(md),
                              delegate = pool[indices(md)])
  
  ## Return MatchedGInteractions object
  return(obj)
}

#' @rdname matchRanges
#' @export
setMethod("matchRanges",
          signature = signature(focal = "GInteractions",
                                pool  = "GInteractions",
                                covar = "formula"),
          definition = matchRanges.MatchedGInteractions)


## Define accessor functions for matchGRanges class --------------------------------------

#' Accessor methods for matchRanges Class
#'
#' @description 
#' Functions that get data from "matchRanges" subclasses
#' such as MatchedDataFrame, MatchedGRanges,
#' and MatchedGInteractions. 
#'
#' @param x Matched object
#' @param ... additional arguments
#'
#' @rdname matchRanges
#' @export
setMethod("focal", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@focal
})

#' @rdname matchRanges
#' @export
setMethod("pool", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool
})

#' @rdname matchRanges
#' @export
setMethod("matched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@delegate
})

#' @rdname matchRanges
#' @export
setMethod("unmatched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool[indices(x, group = "unmatched"),]
})
