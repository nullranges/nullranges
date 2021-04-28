## Generics for Matched class ------------------------------------------------------------

#' Overview
#' 
#' @rdname overview
#' @export
setGeneric("overview", function(x) standardGeneric("overview"))

#' matchedData
#' 
#' @rdname matchedData
#' @export
setGeneric("matchedData", function(x, ...) standardGeneric("matchedData"))

#' Covariates
#' 
#' @rdname covariates
#' @export
setGeneric("covariates", function(x, ...) standardGeneric("covariates"))

#' Indices
#' 
#' @rdname indices
#' @export
setGeneric("indices", function(x, ...) standardGeneric("indices"))

#' Plot
#' 
#' @rdname plot
#' @export
setGeneric("plot")

#' plotCovariates
#' 
#' @rdname plotCovariates
#' @export
setGeneric("plotCovariates", function(x, covar, sets, type, logTransform, ...)
  standardGeneric("plotCovariates"))


## Generics for utils --------------------------------------------------------------------

#' @rdname combnCov
#' @export
setGeneric('combnCov', function(x) standardGeneric('combnCov'))


## Generics for matchedDataFrame/matchedGRanges/matchedGInteractions ---------------------

#' @rdname matchRanges
#' @export
setGeneric("matchRanges", function(focal, pool, covar, method, replace, ...)
  standardGeneric("matchRanges"))

#' Focal
#' 
#' @rdname focal
#' @export
setGeneric("focal", function(x, ...) standardGeneric("focal"))

#' Pool
#' 
#' @rdname pool
#' @export
setGeneric("pool", function(x, ...) standardGeneric("pool"))

#' Matched
#' 
#' @rdname matched
#' @export
setGeneric("matched", function(x, ...) standardGeneric("matched"))

#' Unmatched
#' 
#' @rdname unmatched
#' @export
setGeneric("unmatched", function(x, ...) standardGeneric("unmatched"))
