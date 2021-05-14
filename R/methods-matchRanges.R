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

## Helper function for nearest neighbor matching with replacement
nnMatch <- function(fps, pps, replace) {
  
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

  # TODO we have NOTE here during pkg check: 'no visible binding for global variable'
  
  ## Randomly subsample with replacement among duplicates
  mdt <- amdt[, .(ppsIndex = sample.vec(ppsIndex, N, replace = replace)), by = fps]
  
  ## Add fpsIndex to mdt
  stopifnot(mdt$fps == fpsMap$fps)
  mdt$fpsIndex <- fpsMap$fpsIndex
  
  ## Reorder by fpsIndex (to match fps)
  mdt <- mdt[order(fpsIndex)]
  
  ## Return matched data table
  return(mdt)
  
}

## Helper function to perform rejection sampling
rejectSample <- function(fps, pps) {
  
  ## Ensure fps <= pps
  if (length(fps) > length(pps)) {
    stop("focal must be <= pool for method = 'rejection'.")
  }
  
  ## Kernal density estimates for focal and pool
  df <- kde(fps)
  dg <- kde(pps)
  
  ## Set scale by finding the highest point of density ratios (focal/pool)
  ## This ensures that pool covers focal at all points
  grid <- seq(from=quantile(pps, .001), to=quantile(pps, .999), length=1000)
  scale <- max(predict(df, x=grid) / predict(dg, x=grid))
  
  ## Calculate the probability of accepting each pool
  thresh <- function(x) ifelse(x > 1e-3, x, 0)
  accept_prob <- pmin(thresh(predict(df, x=pps))/(scale * predict(dg, x=pps)), 1)
  
  ## Randomly select pps according to accept probability
  accept <- rbinom(length(pps), size=1, prob=accept_prob)
  
  ## Return random selection equal to the length of pps
  return(accept)
  
}

## Helper function for rejection sampling without replacement
## Not viable for discrete or nearly discrete covariates
rsMatch <- function(fps, pps, replace) {
  
  ## Detect descrete distribution
  accept <- tryCatch(expr = {
    
    ## Rejection sampling
    rejectSample(fps, pps)
    
  }, error = function(e) {
    
    ## Introduce random noise to convert discrete into continuous
    stdev <- seq(1e-06, 0.1, length.out = 100)
    
    ## Iterate over possible noise levels
    for(i in 1:length(stdev)) {
      
      ## Add noise to fps and pps
      nfps <- fps + (rnorm(length(fps), mean = 0, sd = stdev[i]))
      npps <- pps + (rnorm(length(pps), mean = 0, sd = stdev[i]))
      
      ## Bound by 0 and 1
      nfps <- pmax(0, pmin(nfps, 1))
      npps <- pmax(0, pmin(npps, 1))
      
      ## Try rejection sampling
      accept <- suppressWarnings(rejectSample(nfps, npps))
      
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

## Helper function to assign fps and pps to n bins
stratify <- function(fm, pm, n) {
  
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

## Helper function for stratified sampling
ssMatch <- function(fps, pps, replace) {
  
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
    clear = F, total = length(fps) + 1)
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
    pb$update(tokens=list(step=sprintf('Iteration %s, %s bins(s)', i, n)),
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
  pb$update(tokens=list(step=sprintf('Iteration %s, %s bins(s), done!', i, n)),
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

## Helper function that calculates propensity scores
## and implements nearest neighbor matching
propensityMatch <- function(covarData, covars, method, replace) {
  
  ## Assemble covariate formula
  f <- as.formula(paste("id ~", paste(covars, collapse = "+")))
  
  ## Run glm model
  model <- speedglm(formula = f, data = covarData,
                    family = binomial("logit"), fitted = TRUE, model = TRUE)
  
  ## Get propensity scores of focal and pool groups as vectors
  psData <- data.table(ps = predict(model, type = "response"), id = model$model$id)
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
matchRanges.Core <- function(focal, pool, covar, method, replace) {
  
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

## MatchedDataFrame method for matchRanges
matchRanges.MatchedDataFrame <- function(focal, pool, covar, method, replace) {
  
  ## Convert focal and pool to DataFrames
  f <- DataFrame(focal)
  p <- DataFrame(pool)
  
  ## Execute shared GRanges/GInteractions matching code
  md <- matchRanges.Core(f, p, covar, method, replace)
  
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
#' A function that generates a covariate-matched control
#' set of the same datatype.
#'
#' @param focal a DataFrame, GRanges, or GInteractions object containing
#' the focal data to match.
#' @param pool a DataFrame, GRanges, or GInteractions object containing
#' the pool from which to select matches.
#' @param covar a rhs formula with covariates on which to match.
#' @param method character describing which matching method to use.
#' supported options are either 'nearest', 'rejection', or 'stratified'.
#' @param replace logical describing whether to select matches with or without
#' replacement.
#' @param ... additional arguments
#' 
#' @return a covariate-matched control set of data
#'
#' @name matchRanges
#' @rdname matchRanges
#' 
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#'
#' @importFrom rlang f_lhs f_rhs
#' @importFrom speedglm speedglm
#' @importFrom ks kde
#' @import S4Vectors
NULL

#' @rdname matchRanges
#' @export
setMethod("matchRanges",
          signature = signature(focal   = "DF_OR_df_OR_dt",
                                pool    = "DF_OR_df_OR_dt",
                                covar   = "formula",
                                method  = "character",
                                replace = "logical"),
          definition = matchRanges.MatchedDataFrame)

## MatchedGRanges method for matchRanges
matchRanges.MatchedGRanges <- function(focal, pool, covar, method, replace) {
  
  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)
  
  ## Execute matching code to get a Matched object
  md <- matchRanges.Core(f, p, covar, method, replace)
  
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
          signature = signature(focal   = "GRanges",
                                pool    = "GRanges",
                                covar   = "formula",
                                method  = 'character',
                                replace = 'logical'),
          definition = matchRanges.MatchedGRanges)

## MatchedGInteractions method for matchRanges
matchRanges.MatchedGInteractions <- function(focal, pool, covar, method, replace) {
  
  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)
  
  ## Execute matching code to get a Matched object
  md <- matchRanges.Core(f, p, covar, method, replace)
  
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
          signature = signature(focal   = "GInteractions",
                                pool    = "GInteractions",
                                covar   = "formula",
                                method  = 'character',
                                replace = 'logical'),
          definition = matchRanges.MatchedGInteractions)

## Define accessor functions for matchGRanges class --------------------------------------

#' Accessor methods for matchRanges Class
#'
#' Functions that get data from "matchRanges" subclasses
#' such as MatchedDataFrame, MatchedGRanges,
#' and MatchedGInteractions. 
#'
#' @param x Matched object
#' @param ... additional arguments
#' 
#' @rdname matchRangesAccessors
#' @export
setMethod("focal", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@focal
})

#' @rdname matchRangesAccessors
#' @export
setMethod("pool", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool
})

#' @rdname matchRangesAccessors
#' @export
setMethod("matched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@delegate
})

#' @rdname matchRangesAccessors
#' @export
setMethod("unmatched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool[indices(x, group = "unmatched"),]
})
