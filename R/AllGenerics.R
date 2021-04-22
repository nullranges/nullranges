## Generics for Matched class ------------------------------------------------------------

#' @rdname overview
#' @export
setGeneric("overview", function(x) standardGeneric("overview"))

#' @rdname matchedData
#' @export
setGeneric("matchedData", function(x, ...) standardGeneric("matchedData"))

#' @rdname covariates
#' @export
setGeneric("covariates", function(x, ...) standardGeneric("covariates"))

#' @rdname indices
#' @export
setGeneric("indices", function(x, ...) standardGeneric("indices"))

#' @rdname plot
#' @export
setGeneric("plot")

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

#' @rdname focal
#' @export
setGeneric("focal", function(x, ...) standardGeneric("focal"))

#' @rdname pool
#' @export
setGeneric("pool", function(x, ...) standardGeneric("pool"))

#' @rdname matched
#' @export
setGeneric("matched", function(x, ...) standardGeneric("matched"))

#' @rdname unmatched
#' @export
setGeneric("unmatched", function(x, ...) standardGeneric("unmatched"))