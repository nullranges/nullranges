#' @slot matchedData A `data.table` with matched data
#' @slot matchedIndex An integer vector corresponding
#'   to the indices in the `pool` which comprise the
#'   `matched` set.
#' @slot covar A character vector describing the covariates
#'   used for matching.
#' @slot method Character describing replacement method
#'   used for matching.
#' @slot replace TRUE/FALSE describing if matching was
#'   done with or without replacement.
