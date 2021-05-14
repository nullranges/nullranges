## Generics for Matched class ------------------------------------------------------------

#' Overview
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname overview
#' @export
setGeneric("overview", function(x, ...) standardGeneric("overview"))

#' matchedData
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname matchedData
#' @export
setGeneric("matchedData", function(x, ...) standardGeneric("matchedData"))

#' Covariates
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname covariates
#' @export
setGeneric("covariates", function(x, ...) standardGeneric("covariates"))

#' Indices
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname indices
#' @export
setGeneric("indices", function(x, ...) standardGeneric("indices"))

#' @rdname matched-plotting
#' @export
setGeneric("plotCovariates", function(x, covar, sets, type, logTransform, ...)
  standardGeneric("plotCovariates"))

## Generics for utils --------------------------------------------------------------------

#' @param x an object
#' 
#' @rdname combnCov
#' @export
setGeneric("combnCov", function(x) standardGeneric("combnCov"))

## Generics for matchedDataFrame/matchedGRanges/matchedGInteractions ---------------------

#' @rdname matchRanges
#' @export
setGeneric("matchRanges", function(focal, pool, covar, method, replace, ...)
  standardGeneric("matchRanges"))

#' Focal
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname focal
#' @export
setGeneric("focal", function(x, ...) standardGeneric("focal"))

#' Pool
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname pool
#' @export
setGeneric("pool", function(x, ...) standardGeneric("pool"))

#' Matched
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname matched
#' @export
setGeneric("matched", function(x, ...) standardGeneric("matched"))

#' Unmatched
#'
#' @param x an object
#' @param ... additional arguments
#' 
#' @rdname unmatched
#' @export
setGeneric("unmatched", function(x, ...) standardGeneric("unmatched"))
