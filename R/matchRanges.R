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
matchRanges <- function(x, univ, covar, method = "nearest") {

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

# matchRangesOld2 <- function(x, univ, covar, method = "nearest") {
#   
#   ## Ensure covar columns are shared between x and univ
#   if(!(all(covar %in% colnames(mcols(x))) & 
#        all(covar %in% colnames(mcols(univ))))) {
#     stop("covar must be columns in both x and univ.", call. = T)
#   }
#   
#   ## Create data frame for matchit
#   df <- as.data.frame(cbind(id = factor(c(rep(1, length(x)), rep(0, length(univ)))),
#                             rbind(mcols(x)[covar], mcols(univ)[covar])))
#   
#   ## Assemble covariate formula
#   f <- as.formula(paste("id ~", paste(covar, collapse = "+")))
#   
#   ## Create matched control object
#   matchObj <- MatchIt::matchit(data = df, formula = f, method = method)
#   
#   ## Extract matched indicies
#   matchIndex <- as.numeric(as.vector(matchObj$match.matrix))
#   
#   ## Return matchced controls
#   return(c(x, univ)[matchIndex])
#   
# }

# matchRangesOld <- function(x, univ, covar, method = "nearest") {
# 
#   ## Add identifier for focal & universe groups
#   mcols(x)$mrid <- 1
#   mcols(univ)$mrid <- 0
# 
#   ## Combine & convert to data.frame
#   df <- as.data.frame(c(x, univ))
# 
#   ## Ensure covar columns are shared between x and univ
#   if(!(all(covar %in% colnames(mcols(x))) &
#        all(covar %in% colnames(mcols(univ))))) {
#     stop("covar must be columns in both x and univ.", call. = T)
#   }
# 
#   ## Assemble covariate formula
#   f <- as.formula(paste("mrid ~", paste(covar, collapse = "+")))
# 
#   ## Create matched control object
#   matchObj <- MatchIt::matchit(data = df, formula = f, method = method)
# 
#   ## Extract matched data
#   md <- MatchIt::match.data(matchObj)
# 
#   ## Convert to GRanges, filter for controls, remove Matchit columns, restore seqinfo
#   ctls <-
#     md %>%
#     plyranges::as_granges(keep_mcols = T) %>%
#     plyranges::filter(mrid == 0) %>%
#     plyranges::select(-c("mrid", "distance", "weights", "subclass")) %>%
#     `seqinfo<-`(value = seqinfo(univ))
# 
#   ## Return matched controls
#   return(ctls)
# 
# }

