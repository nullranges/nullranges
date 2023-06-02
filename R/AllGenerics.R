## Generics for Matched class ------------------------------------------------------------

#' @rdname overview
#' @export
setGeneric("overview", function(x, digits = 2, ...) standardGeneric("overview"))

#' @rdname matchedData
#' @export
setGeneric("matchedData", function(x, ...) standardGeneric("matchedData"))

#' @rdname covariates
#' @export
setGeneric("covariates", function(x, ...) standardGeneric("covariates"))

#' @rdname method
#' @export
setGeneric("method", function(x, ...) standardGeneric("method"))

#' @rdname withReplacement
#' @export
setGeneric("withReplacement", function(x, ...) standardGeneric("withReplacement"))

#' @rdname indices
#' @export
setGeneric("indices", function(x, set = 'matched', ...) standardGeneric("indices"))

#' @rdname plotPropensity
#' @export
setGeneric("plotPropensity", function(x, 
                                      sets = c('focal',
                                               'matched',
                                               'pool',
                                               'unmatched'),
                                      type = NULL,
                                      log = NULL,
                                      ...)
  standardGeneric("plotPropensity"))

#' @rdname plotCovariate
#' @export
setGeneric("plotCovariate", function(x,
                                     covar = NULL,
                                     sets = c('focal',
                                              'matched',
                                              'pool',
                                              'unmatched'),
                                     type = NULL,
                                     log = NULL,
                                     ...)
  standardGeneric("plotCovariate"))

## Generics for utils --------------------------------------------------------------------

#' @rdname combnCov
#' @export
setGeneric("combnCov", function(x, ...) standardGeneric("combnCov"))

#' @rdname makeExampleMatchedDataSet
#' @export
setGeneric("makeExampleMatchedDataSet", function(type = 'DataFrame',
                                                 matched = FALSE,
                                                 method = 'rejection', 
                                                 replace = FALSE,
                                                 ...)
  standardGeneric("makeExampleMatchedDataSet"))

## Generics for matchedDataFrame/matchedGRanges/matchedGInteractions ---------------------

#' @rdname matchRanges
#' @export
setGeneric("matchRanges", function(focal, pool, covar,
                                   method = 'nearest',
                                   replace = TRUE,
                                   ...)
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
