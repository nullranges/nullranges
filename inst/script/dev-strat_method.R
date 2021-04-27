## Load packages and data ----------------------------------------------------------------

library(nullranges)
library(nullrangesData)
library(hictoolsr)
library(magrittr)
library(GenomicRanges)
library(speedglm)


## Load enh-prom pairs (and shorten the name)
data("enhPromContactFreqHg19")
epp <- enhPromContactFreqHg19

## Load loops
loops <- fread(system.file("extdata/hic/MOMA_SIP_10kbLoops_Merged.txt",
                           package = 'nullrangesData')) %>% as_ginteractions()

## Annotate looped and unlooped enh-prom pairs
epp$loopedEP <- FALSE
epp$loopedEP[countOverlaps(epp, loops) > 0] <- TRUE

## Define focal, pool and covar
focal <- mcols(epp[epp$epDistance <= 40e03])
pool  <- mcols(epp[epp$epDistance >= 40e03])
covar <- ~loopedEP
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

## Develop new stratified matching method

