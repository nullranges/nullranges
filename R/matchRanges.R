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
#' @param method Character string passed to MatchIt used to calculate matches
#' 
#' @return Returns a GRanges object of matched controls
#' 
#' @examples 
#' 
#' matchRanges(x = up_ov, univ = ctl_ov, covar = c("peakStrength", "gc"))
#' 
#'
#' @export
#' 
matchRanges <- function(x, univ, covar, method = "nearest") {
  
  ## Add identifier for focal & universe groups
  mcols(x)$mrid <- 1 
  mcols(univ)$mrid <- 0 
  
  ## Combine & convert to data.frame
  df <- as.data.frame(c(x, univ))
  
  ## Ensure covar columns are shared between x and univ
  if(!(all(covar %in% colnames(mcols(x))) & 
       all(covar %in% colnames(mcols(univ))))) {
    stop("covar must be columns in both x and univ.", call. = T)
  }
  
  ## Assemble covariate formula
  f <- as.formula(paste("mrid ~", paste(covar, collapse = "+")))
  
  ## Create matched control object
  matchObj <- MatchIt::matchit(data = df, formula = f, method = method)
  
  ## Extract matched data
  md <- match.data(matchObj)
  
  ## Convert to GRanges, filter for controls, remove Matchit columns, restore seqinfo
  ctls <- 
    md %>%
    plyranges::as_granges(keep_mcols = T) %>%
    plyranges::filter(mrid == 0) %>%
    plyranges::select(-c("mrid", "distance", "weights", "subclass")) %>%
    `seqinfo<-`(value = seqinfo(univ))
  
  ## Return matched controls
  return(ctls)
  
}

