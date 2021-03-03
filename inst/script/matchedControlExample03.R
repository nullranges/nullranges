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

## Annotate looped enh-prom pairs
loopedEP <- subsetByOverlaps(enhPromPairs, loopGR)
unLoopedEP <- subsetByOverlaps(enhPromPairs, loopGR, invert = T)

## Generate matched controls
system.time({
  matched <- matchRanges(x = loopedEP,
                         univ = unLoopedEP,
                         covar = c("epDistance", "anchor1.peakStrength"))  
})


## See how well the distributions are matched
plot(density(loopedEP$epDistance))
lines(density(matched$epDistance), col = "blue", lty=2, lwd = 3)
lines(density(s$epDistance), col = "green", lty=2, lwd = 3)

plot(density(log2(loopedEP$anchor1.peakStrength)), ylim = c(0, 0.4))
lines(density(log2(matched$anchor1.peakStrength)), col = "blue", lty=2, lwd = 3)
lines(density(log2(s$anchor1.peakStrength)), col = "green", lty=2, lwd = 3)

## Visualize the correlation between looped and unlooped enh-prom
barplot(c(UnLooped = cor(s$anchor1.peakFC, s$anchor2.log2FoldChange),
          MatchedControl = cor(matched$anchor1.peakFC, matched$anchor2.log2FoldChange),
          Looped = cor(loopedEP$anchor1.peakFC, loopedEP$anchor2.log2FoldChange)),
        ylab = "Correlation",
        main = "Correlation between FC(H3K27Ac) and FC(RNA)")


## Differential enhancer at one end and boxplot of RNA at the other end

# library(ggplot2)
# ggplot(enhPromDt, aes(x = first.peakFC, y = second.log2FoldChange))+
#   facet_grid(~looped)+
#   geom_hex()