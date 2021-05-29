## Exported general utility functions ----------------------------------------------------

## Function for creating combinations of covariates
combnCovariates <- function(x) {

  ## Generate covariate list for expanding
  covarList <- lapply(x, function(y) c(y, 0, 0))

  ## Generate all unique covariate combinations
  uniqCovComb <- data.table(unique(expand.grid(covarList)))

  ## Combine by row into a formula
  f <- paste0("~", apply(uniqCovComb, 1, paste0, collapse = '+'))

  ## Remove 0's and the last formula (always empty)
  f <- gsub("\\+0|0\\+", "", f[-length(f)])

  ## Return formulae
  return(f)

}

#' Function for creating combinations of covariates
#'
#' @param x Character vector of covariates to combine
#'
#' @return Returns a character vector of formulae combinations.
#'
#' @examples
#' combnCov(x = c('a', 'b', 'c'))
#'
#' @rdname combnCov
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @export
setMethod("combnCov", signature(x="character"), combnCovariates)

## Function for generating an example dataset
make_example_matched_data_set <- function(type, matched, method, replace) {
  
  ## Parse type argument
  type <- match.arg(type, choices = c('data.frame',
                                      'data.table',
                                      'DataFrame',
                                      'GRanges',
                                      'GInteractions'))
  
  ## Generate example covariate data
  set.seed(123)
  df <- data.frame(
    treated = c(rep(TRUE, 500),
                rep(FALSE, 1e4)),
    covar1 = c(abs(rnorm(500, mean = 4, sd = 2)),
               runif(1e4, min = 0, max = 12)),
    covar2 = c(sample(letters[seq_len(5)],
                      size = 500,
                      replace = TRUE,
                      prob = c(0.1, 0.3, 0.4, 0.1, 0.05)),
               sample(letters[seq_len(5)],
                      size = 1e4,
                      replace = TRUE,
                      prob = c(0.4, 0.3, 0.1, 0.1, 0.05)))
  )
  df$covar2 <- as.factor(df$covar2) # remove after fixing overview()

  ## Generate example data.frame/data.table/DataFrame
  if (identical(type, 'data.frame')) out <- df
  if (identical(type, 'data.table')) out <- as.data.table(df)
  if (identical(type, 'DataFrame')) out <- DataFrame(df)
  
  ## Generate example GRanges
  if (identical(type, 'GRanges')) {
    out <- GRanges(seqnames = 'chr1', 
                   ranges = IRanges(start = seq_len(nrow(df)),
                                    width = 100))
    
    mcols(out) <- df
  }
  
  ## Generate example GInteractions
  if (identical(type, 'GInteractions')) {
    gr <- GRanges(seqnames = 'chr1', 
                  ranges = IRanges(start = seq_len(nrow(df)),
                                   width = 100))
    
    out <- GInteractions(anchor1 = seq_len(nrow(df)),
                         anchor2 = seq_len(nrow(df)),
                         regions = gr)
    
    mcols(out) <- df
  }
  
  ## Return dataset matched or not
  if (matched) {
    out <- matchRanges(focal = out[out$treated,],
                       pool = out[!out$treated,],
                       covar = ~covar1 + covar2,
                       method = method,
                       replace = replace)
  }
  
  return(out)
}


#' Function for generating an example matchRanges or Matched dataset
#' 
#' This function will generate an example dataset as either 1) input
#' for `matchRanges()` (when `matched = TRUE`) or 2) a 
#' Matched Object (when `matched = FALSE`).
#'
#' @param type Character designating which type of dataset to make.
#'   options are one of 'data.frame', 'data.table', 'DataFrame',
#'   'GRanges', or 'GInteractions'.
#' @param matched TRUE/FALSE designating whether to generate a
#'   Matched dataset (`matched = TRUE`) or an input dataset
#'   for `matchRanges()` (`matched = FALSE`).
#' @param method Character describing which matching method to use.
#'   Supported options are either 'nearest', 'rejection', or 'stratified'.
#' @param replace TRUE/FALSE describing whether to select matches with
#'   or without replacement.
#' @param ... Additional arguments
#'
#' @return Returns an example Matched dataset or an example dataset for
#'   input to `matchRanges()`.
#'
#' @examples
#' ## Make examples for matchRanges() (i.e matched = FALSE)
#' makeExampleMatchedDataSet()
#' head(makeExampleMatchedDataSet(type = 'data.frame', matched = FALSE))
#' makeExampleMatchedDataSet(type = 'data.table', matched = FALSE)
#' makeExampleMatchedDataSet(type = 'DataFrame', matched = FALSE)
#' makeExampleMatchedDataSet(type = 'GRanges', matched = FALSE)
#' makeExampleMatchedDataSet(type = 'GInteractions', matched = FALSE)
#'
#' ## Make Matched class examples (i.e. matched = TRUE)
#' makeExampleMatchedDataSet(matched = TRUE)
#' makeExampleMatchedDataSet(type = 'DataFrame', matched = TRUE,
#'                           method = 'rejection',
#'                           replace = FALSE)
#' makeExampleMatchedDataSet(type = 'GRanges', matched = TRUE,
#'                           method = 'rejection',
#'                           replace = FALSE)
#' makeExampleMatchedDataSet(type = 'GInteractions', matched = TRUE,
#'                           method = 'rejection',
#'                           replace = FALSE)
#'
#' @rdname makeExampleMatchedDataSet
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @import GenomicRanges
#' @import InteractionSet
#' @import S4Vectors
#' @export
setMethod("makeExampleMatchedDataSet",
          signature = signature(type = 'character_OR_missing',
                                matched = 'logical_OR_missing',
                                method = 'character_OR_missing',
                                replace = 'logical_OR_missing'),
          definition = make_example_matched_data_set)