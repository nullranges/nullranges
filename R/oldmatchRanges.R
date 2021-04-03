#' Generate a set of matched control GRanges
#'
#' @description
#' Given an input of focal GRanges, universe GRanges and covariates, this function will
#' generate a set of GRanges selected from the universe that are matched on the covariates
#' to the focal ranges. This function wraps the covariate matching function from the
#' package 'MatchIt'.
#'
#' @param x GRanges object containing focal ranges to match
#' @param univ GRanges object containing the universe of ranges from which to select
#' @param covar Character vector of covariates to match on
#' 
#' @return Returns a GRanges object of matched controls
#'
#' @examples
#'
#' matchRanges(x = up_ov, univ = ctl_ov, covar = c("peakStrength", "gc"))
#' @export
#'
oldmatchRanges <- function(x, univ, covar, method = "nearest") {

  ## Create data frame for matching covariates
  df <- as.data.frame(cbind(id = factor(c(rep(1, length(x)), rep(0, length(univ)))),
                            rbind(mcols(x)[covar], mcols(univ)[covar])))
  
  ## Assemble covariate formula
  f <- as.formula(paste("id ~", paste(covar, collapse = "+")))
  
  ## Run glm model
  model <- glm(formula = f, data = df, family = binomial("logit"))
  
  ## Get x and univ propensity scores as vectors
  ps_df <- data.frame(ps = predict(model, type = "link"), id = model$model$id)
  xps <- ps_df[ps_df$id == 1,1]
  ups <- ps_df[ps_df$id == 0,1]
  
  ## Create data table with original ups index and setkey to sort
  dt <- data.table(ups, val = ups, index = 1:length(ups))
  setkey(dt, ups)
  
  ## Find the ups index of the nearest neighbor to each value in xps with a rolling join
  nnI <- dt[.(xps), index, roll = 'nearest', mult = 'first', with = T]
  
  ## Return the matched controls
  return(univ[nnI])

}
