#' @section Accessor methods for Matched Class:
#' Functions that get data from Matched subclasses (`x`)
#' such as MatchedDataFrame, MatchedGRanges,
#' and MatchedGInteractions are listed below:
#' \itemize{
#'   \item{`matchedData(x)`: }{Get matched data from a Matched object}
#'   \item{`covariates(x)`: }{Get covariates from a Matched object}
#'   \item{`method(x)`: }{Get matching method used for Matched object}
#'   \item{`withReplacement(x)`: }{Get replace method}
#'   \item{`indices(x, set)`: }{Get indices of matched set}
#' }
#' For more detail check the help pages for these functions.
#' 
#' @examples 
#' ## Make MatchedDataFrame example
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' 
#' ## Accessor functions for Matched class
#' matchedData(x)
#' covariates(x)
#' method(x)
#' withReplacement(x)
#' head(indices(x, set = 'matched'))
#' 
#' @seealso [matchedData], [covariates], [method],
#'   [withReplacement], [indices]
