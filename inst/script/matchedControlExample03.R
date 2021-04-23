## More relevant motivating example:
##
## Is fold-change in k27ac (macrophage/monocyte or PMA/NON) at enhancers correlated with
## RNA when enh-pro is looped vs. not looped?
##
## How about when matched for the distance of the enhancer to its gene promoter or the
## K27ac signal of the enhancer?

## Loading libraries and data ------------------------------------------------------------

## Load required libraries
library(data.table)
library(nullranges)
library(nullrangesData)
library(plyranges)
library(InteractionSet)
library(strawr)


## Load RNA-seq differential genes & H3K27Ac peak counts
data("monoMacroDEGHg19")
data("H3K27acPeakCountsHg19")

## Assemble enhancer-promoter pairs ------------------------------------------------------

## Define gene promoters
prom <- promoters(monoMacroDEGHg19, upstream = 2000, downstream = 200)

## Define enhancers (remove K27 signal that falls within a promoter)
enh <- subsetByOverlaps(H3K27acPeakCountsHg19, prom, invert = T)

## Define 2Mb regions around gene promoters
promRegions <- 
  prom %>%
  resize(., width = 2e06, fix = "center") %>%
  suppressWarnings() %>%
  trim()

## Generate enhancer-promoter pairs for enhancers that fall within 2Mb of a gene promoter
enhPromPairs <-
  findOverlaps(enh, promRegions, type = "within") %>%
  Pairs(enh, prom, hits = .)

## Convert to GInteractions object (while adding mcols)
enhPromPairs <- GInteractions(first(enhPromPairs), second(enhPromPairs))

## Calculate the distance between enh-pro pairs
enhPromPairs$epDistance <- pairdist(enhPromPairs, type="gap")


## Reassign enh-prom pairs to hic res bins -----------------------------------------------

## Extract enh and prom
enhAnchor <- first(enhPromPairs)
proAnchor <- second(enhPromPairs)

## Define resolution to bin by
res <- 10000

## Resize by the center of the enhancer and bin
enhAnchor %<>%
  resize(width = 0, fix = 'center') %>%
  mutate(start = floor(start/res)*res,
         end = floor(start/res)*res + res - 1)

## Resize at the TSS of promoter and bin
proAnchor %<>%
  resize(width = 0, fix = 'end') %>%
  mutate(start = start - 200, end = start) %>%
  mutate(start = floor(start/res)*res,
         end = floor(start/res)*res + res - 1)

## Combine into binned enhancer-promoter pairs
binnedEP <- GInteractions(enhAnchor, proAnchor)

## Add back metadata columns
mcols(binnedEP) <- mcols(enhPromPairs)


## Collect interaction frequencies between enhancer-promoter pairs -----------------------

## TODO
## Write helper function to bin bedpe interactions by resolution
##    Make it flexible enough to handle data.tables/frames and GInteraction objects
## Apply helper functino in extract Counts
## Move this function to JuicePouch repository

#' @param bedpe GInteraction object that has been binned to the correct resolution
#' @param hic Path to .hic file
#' @param res Resolution of bedpe bins
#' @param chroms Vector of chromosomes to extract
#' 
extractCounts <- function(bedpe, hic, chroms = c(1:22, 'X', 'Y'), res = 10000) {

  ## Start progress bar
  pb <- progress::progress_bar$new(
    format = "  :step [:bar] :percent elapsed: :elapsedfull",
    clear = F, total = length(hic)*length(chroms)*3+1)
  pb$tick(0)

  ## Convert chroms to character
  chroms <- as.character(chroms)
  
  ## Split binned binnedEP by chromosome
  chrBedpe <- lapply(chroms, function(chr) {
    bedpe[seqnames(first(bedpe)) == paste0('chr',chr)]
  })
  
  ## Add names
  names(chrBedpe) <- chroms
  
  ## Loop through each hic file and chromosome
  for(i in 1:length(hic)) {
    for(j in 1:length(chroms)) {
      
      ## Update progress
      pb$tick(tokens=list(step=sprintf('Extracting chr%s from %s',
                                       chroms[j], basename(hic[i]))))
      
      ## Construct straw format from binnedEP
      ep <- data.table(
        x = start(first(chrBedpe[[chroms[j]]])),
        y = start(second(chrBedpe[[chroms[j]]])),
        counts = 0
      )
      
      ## Swap row assignment for out of order x and y
      ep[y < x, `:=`(y,x)]
      
      ## Pull out the entire chromosome sparse matrix
      sparseMat <- as.data.table(strawr::straw(norm = "NONE", fname = hic[i],
                                               chr1loc = chroms[j],
                                               chr2loc = chroms[j],
                                               unit = "BP",
                                               binsize = res,
                                               matrix = "observed"))
      
      ## Update progress
      pb$tick(tokens=list(step=sprintf('Extracting chr%s from %s',
                                       chroms[j], basename(hic[i]))))
      
      ## Set keys
      setkeyv(sparseMat, c('x', 'y'))
      
      ## Get counts by key
      ep$counts <- sparseMat[ep]$counts
      
      ## Set unmatched counts (NA) to 0
      ep[is.na(counts), counts := 0]
      
      ## Add counts back to chrBedpe
      mcols(chrBedpe[[chroms[j]]])[basename(hic[i])] <- ep$counts
      
      ## Update progress
      pb$tick(tokens=list(step=sprintf('Extracting chr%s from %s',
                                       chroms[j], basename(hic[i]))))
      
    }
  }
  
  ## Combine results
  chrBedpe <- unname(chrBedpe)
  combinedResults <- do.call(c, chrBedpe)
  
  ## Close progress bar
  pb$tick(tokens = list(step = "Done!"))
  if(pb$finished) pb$terminate()
  
  ## Return combined results
  return(combinedResults)
}


## Extract loop counts for each biorep hic file
result <- extractCounts(bedpe = binnedEP,
                        hic = c("../hic/CI_THP1_O_1_inter_30.hic",
                                "../hic/CI_THP1_O_2_inter_30.hic",
                                "../hic/CI_THP1_A_1_inter_30.hic",
                                "../hic/CI_THP1_A_2_inter_30.hic"),
                        chroms = c(1:22,'X','Y'),
                        res = 10000)

result



## Load loops ----------------------------------------------------------------------------

## Load loop BEDPE
loops <- fread(system.file("extdata/hic/MOMA_SIP_10kbLoops_Merged.txt",
                           package = 'nullrangesData'))

## Define anchors
a1 <- loops[,1:3] %>% `colnames<-`(c("seqnames", "start", "end")) %>% as_granges()
a2 <- loops[,4:6] %>% `colnames<-`(c("seqnames", "start", "end")) %>% as_granges()

## Create GInteractions object
loopGR <- GInteractions(anchor1 = a1, anchor2 = a2)


## Perform matched control comparisons ---------------------------------------------------

## Add contact frequency to enhPromPairs
enhPromPairs$contactFreq <- mcols(result)[grep("hic", colnames(mcols(result)))] %>%
  as.matrix() %>%
  rowSums()


## Annotate looped enh-prom pairs
loopedEP <- subsetByOverlaps(enhPromPairs, loopGR)
unLoopedEP <- subsetByOverlaps(enhPromPairs, loopGR, invert = T)


## Generate matched controls
system.time({
  matched <- matchRanges(x = loopedEP,
                         univ = unLoopedEP,
                         covar = ~epDistance+anchor1.peakStrength+contactFreq,
                         method = 'nearest', replace = TRUE)  
})

## Diagnostic plots: see how well the distributions are matched
plot(density(loopedEP$epDistance))
lines(density(matched$epDistance), col = "blue", lty=2, lwd = 3)
lines(density(unLoopedEP$epDistance), col = "green", lty=2, lwd = 3)

plot(density(log2(loopedEP$anchor1.peakStrength)), ylim = c(0, 0.4))
lines(density(log2(matched$anchor1.peakStrength)), col = "blue", lty=2, lwd = 3)
lines(density(log2(unLoopedEP$anchor1.peakStrength)), col = "green", lty=2, lwd = 3)

plot(density((loopedEP$contactFreq)), ylim = c(0, 0.00065), xlim = c(-2000, 15000))
lines(density((matched$contactFreq)), col = "blue", lty=2, lwd = 3)
lines(density((unLoopedEP$contactFreq)), col = "green", lty=2, lwd = 3)

plot(density(log2(loopedEP$contactFreq+1)), ylim = c(0, 0.4))
lines(density(log2(matched$contactFreq+1)), col = "blue", lty=2, lwd = 3)
lines(density(log2(unLoopedEP$contactFreq+1)), col = "green", lty=2, lwd = 3)

## Visualize trends ----------------------------------------------------------------------

## Visualize the correlation between looped and unlooped enh-prom
barplot(c(UnLooped = cor(unLoopedEP$anchor1.peakFC, unLoopedEP$anchor2.log2FoldChange),
          MatchedControl = cor(matched$anchor1.peakFC, matched$anchor2.log2FoldChange),
          Looped = cor(loopedEP$anchor1.peakFC, loopedEP$anchor2.log2FoldChange)),
        ylab = "Correlation",
        main = "Correlation between FC(H3K27Ac) and FC(RNA)")


## Differential enhancer at one end and boxplot of RNA at the other end
Log2FC_RNA <- list(
  unLooped = unLoopedEP$anchor2.log2FoldChange[unLoopedEP$anchor1.peakFC >= 2],
  matched = matched$anchor2.log2FoldChange[matched$anchor1.peakFC >= 2],
  looped = loopedEP$anchor2.log2FoldChange[loopedEP$anchor1.peakFC >= 2]
)

library(ggplot2)
df <- data.frame(group = factor(c(rep("No E-P Loop", length(Log2FC_RNA$unLooped)),
                                  rep("Matched Control", length(Log2FC_RNA$matched)),
                                  rep("Looped E-P", length(Log2FC_RNA$looped))),
                                levels = c("No E-P Loop", "Matched Control", "Looped E-P")),
                 Log2FC_RNA = c(Log2FC_RNA$unLooped,
                                Log2FC_RNA$matched,
                                Log2FC_RNA$looped))
  
ggplot(df, aes(x = group, y = Log2FC_RNA, color = group))+
  geom_boxplot(outlier.colour = NA) +
  ylim(c(-5, 5))+
  theme_bw()





## Work on visualizations to assess matching quality -------------------------------------

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


plot(jitter(rep(1.5, length(loopedEP$epDistance)))~
       jitter(loopedEP$epDistance), ylim = c(0, 2), cex = 0.1)
points(jitter(rep(1.0, length(matched$epDistance)))~
         jitter(matched$epDistance), cex = 0.1, col = "blue")
points(jitter(rep(0.5, length(unLoopedEP$epDistance)))~
         jitter(unLoopedEP$epDistance), cex = 0.1, col = "green")