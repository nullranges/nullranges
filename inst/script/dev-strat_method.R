## Load packages and data ----------------------------------------------------------------

library(nullranges)
library(nullrangesData)
library(hictoolsr)
library(magrittr)
library(GenomicRanges)
library(speedglm)
library(data.table)


## Load enh-prom pairs (and shorten the name)
data("enhPromContactFreqHg19")
epp <- enhPromContactFreqHg19

## Load loops
loops <- fread(system.file("extdata/hic/MOMA_SIP_10kbLoops_Merged.txt",
                           package = 'nullrangesData')) %>% as_ginteractions()

## Annotate looped and unlooped enh-prom pairs
epp$loopedEP <- FALSE
epp$loopedEP[countOverlaps(epp, loops) > 0] <- TRUE

## Prepare data for new method -----------------------------------------------------------

## Define sample.vec to handle vectors of varying length
sample.vec <- function(x, ...) x[sample(length(x), ...)]

## Define focal, pool and covar
focal <- mcols(epp[epp$epDistance <= 40e03])
pool  <- mcols(epp[epp$epDistance >= 40e03])
# covar <- ~loopedEP
# covar <- ~contactFreq
# covar <- ~anchor1.peakStrength
# covar <- ~loopedEP + contactFreq
covar <- ~loopedEP + contactFreq + anchor1.peakStrength

# focal <- mcols(epp[epp$loopedEP == TRUE])
# pool  <- mcols(epp[epp$loopedEP == FALSE])
# covar <- ~epDistance
# covar <- ~contactFreq
# covar <- ~anchor1.peakStrength
# covar <- ~epDistance + contactFreq
# covar <- ~epDistance + contactFreq + anchor1.peakStrength

method <- 'stratified'
replace <- FALSE

## Extract covariates from formula as character vector
covars <- nullranges:::parseFormula(covar)

## Check that all covariates are in both focal and pool
if (!(all(covars %in% colnames(focal)) &
      all(covars %in% colnames(pool)))) {
  stop("All variables in covar must be columns in both focal and pool.")
}

## Check method and replace arguments
method <- match.arg(method, choices = c('nearest', 'rejection', 'stratified'))

if (isFALSE(replace) & nrow(focal) >= nrow(pool)) 
  stop("focal must be <= pool when replace = FALSE.")

if (method == 'nearest' & isFALSE(replace))
  stop("nearest neighbor matching without replacement not available.")

## Create data table with covariate data
covarData <- as.data.table(cbind(id = factor(c(rep(1, nrow(focal)),
                                               rep(0, nrow(pool)))),
                                 rbind(focal[covars], pool[covars])))

## Assemble covariate formula
f <- as.formula(paste("id ~", paste(covars, collapse = "+")))

## Run glm model
model <- speedglm(formula = f, data = covarData,
                  family = binomial("logit"), fitted = TRUE, model = TRUE)

## Get propensity scores of focal and pool groups as vectors
psData <- data.table(ps = predict(model, type = "response"), id = model$model$id)
fps <- psData[id == 1, ps]
pps <- psData[id == 0, ps]

## Add propensity scores to covarData
covarData$ps <- psData$ps

## Develop new stratified matching method ------------------------------------------------

## Helper function to assign fps and pps to n bins
stratify <- function(fm, pm, n) {
  
  ## Define breaks using fps and pps
  mn <- min(c(fm$fps, pm$pps))
  mx <- max(c(fm$fps, pm$pps))
  br <- seq(from=mn, to=mx, by=(mx-mn)/n)
  
  ## Assign fps and pps to bin
  fm$bin <- .bincode(fm$fps, br, include.lowest = TRUE)
  pm$bin <- .bincode(pm$pps, br, include.lowest = TRUE)
  
  ## Assign indices to bins
  fpsBins <- fm[, .(fpsN = .N, fpsIndices = list(fpsIndex)), by = bin]
  ppsBins <- pm[, .(ppsN = .N, ppsIndices = list(ppsIndex)), by = bin]
  
  ## Define strata by joining fps and pps on bins
  strata <- fpsBins[ppsBins, on = 'bin']
  
  return(strata)
  
}

## Rewrite best-fit method ####

system.time({
  ## Initialize results, fpsOptions and ppsOptions
  results <- data.table(bin=integer(), fpsIndex=integer(), ppsIndex=integer())
  fpsOptions <- data.table(fps, val = fps, fpsIndex = seq_along(fps))
  ppsOptions <- data.table(pps, val = pps, ppsIndex = seq_along(pps))
  
  ## Set flags and iteration count
  skip <- FALSE
  i <- 1
  
  while (nrow(results) != nrow(focal)) {
    
    ## Update n
    if (skip)  n <- floor(n/2)
    if (!skip) n <- length(unique(c(fpsOptions$fps, ppsOptions$pps)))
    
    ## Stratify ps by bins and match focal and pool
    strata <- stratify(fpsOptions, ppsOptions, n)
    
    while (nrow(strata[!is.na(fpsN) & fpsN <= ppsN]) == 0) {
      ## Enter faster n searching
      skip <- TRUE
      n <- floor(n/2)
      
      ## Stratify ps by bins and match focal and pool
      strata <- stratify(fpsOptions, ppsOptions, n)
    }
    
    ## Print out step
    print(sprintf("iteration %s: %s %% complete, %s bin(s)", i,
                  round(nrow(results)/nrow(focal) * 100, 2), n))
    i <- i + 1
    
    ## Assign indices that can be sampled
    set.seed(123)
    result <-
      strata[!is.na(fpsN) & fpsN <= ppsN,
             .(fpsIndex = unlist(fpsIndices),
               ppsIndex = sample.vec(unlist(ppsIndices), fpsN, replace = replace)),
             by = bin]
    
    ## Append to results
    results <- rbind(results, result)
    
    ## Remove assigned indices from options
    fpsOptions <- fpsOptions[!fpsIndex %in% result$fpsIndex]
    ppsOptions <- ppsOptions[!ppsIndex %in% result$ppsIndex]
    
  }
  
  ## Print out step
  print(sprintf("iteration %s: %s %% complete, %s bin(s)", i,
                round(nrow(results)/nrow(focal) * 100, 2), n))
  i <- i + 1
})

## Reorder by fpsIndex
results <- results[order(results$fpsIndex)]

## Assemble matched data table
mdt <- data.table(fps = fps[results$fpsIndex],
                  ppsIndex = results$ppsIndex,
                  fpsIndex = results$fpsIndex)

## Assemble information by group
matchedData <- rbind(
  covarData[id == 1, c(.SD, group = 'focal')],
  covarData[id == 0][mdt$ppsIndex, c(.SD, group = 'matched')],
  covarData[id == 0, c(.SD, group = 'pool')],
  covarData[id == 0][!mdt$ppsIndex, c(.SD, group = 'unmatched')]
)

## Matched indicies
matchedIndex <- mdt$ppsIndex



## Combine information in Matched object
obj <- nullranges:::Matched(matchedData = matchedData,
                            matchedIndex = matchedIndex,
                            covar = covars)

## Look at results -----------------------------------------------------------------------

mean(matchedData(obj)[group == 'focal', ps] -
       matchedData(obj)[group == 'matched', ps])

overview(obj)
covariates(obj)

table(matchedData(obj)[group == 'focal']$loopedEP)
table(matchedData(obj)[group == 'matched']$loopedEP)

length(table(nullranges::indices(obj, 'matched'))[table(nullranges::indices(obj, 'matched'))>1])

plot(obj, type = 'lines')
plot(obj, type = 'lines') + ggplot2::xlim(c(0, 0.02))
plot(obj, type = 'lines') + ggplot2::scale_x_log10()

plotCovariates(obj, type = 'lines')
plotCovariates(obj, type = 'lines', logTransform = TRUE)

plotCovariates(obj, covar = "epDistance", type = 'lines', logTransform = FALSE)
plotCovariates(obj, covar = "contactFreq", type = 'lines', logTransform = TRUE)
plotCovariates(obj, covar = "anchor1.peakStrength", type = 'lines', logTransform = TRUE)
