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
  f <- f[-length(f)] %>% gsub("\\+0|0\\+", "", .)
  
  ## Return formulae
  return(f)
  
}

#' Function for creating combinations of covariates
#' 
#' @param x character vector of covariates to combine
#' 
#' @return returns a character vector of formulae combinations.
#' 
#' @examples 
#' 
#' combnCov(x = c('a', 'b', 'c'))
#' 
#' @rdname combnCov
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @import magrittr
#' @export
setMethod("combnCov", signature(x="character"), combnCovariates)


