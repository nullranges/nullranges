## Example using matching function(s)


## Loading libraries and data ------------------------------------------------------------

## Load required libraries
library(nullranges)
library(nullrangesData)
library(magrittr)
library(plyranges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

## Load RNA-seq differential genes & H3K27Ac peak counts
data("monoMacroDEGHg38")
data("H3K27acPeakCountsHg38")

## Simplify dataset names
de_genes <- monoMacroDEGHg38
chipPeaks <- H3K27acPeakCountsHg38

## Basic processing ----------------------------------------------------------------------

## Get promoters of de_genes
de_gene_prom <- promoters(x = de_genes, upstream = 2000, downstream = 200)

## Define up, down, & control promoters
up_gene_prom  <- de_gene_prom %>% filter(log2FoldChange >= 2 & padj < 0.01)
dn_gene_prom  <- de_gene_prom %>% filter(log2FoldChange <= -2 & padj < 0.01)
ctl_gene_prom <- de_gene_prom %>% filter(!(padj < 0.01))
# ctl_gene_prom <- subsetByOverlaps(x = de_gene_prom,
#                                   ranges = c(up_gene_prom, dn_gene_prom), invert = T)

## Find k27 overlap with up/ctl promoters
up_ov  <- subsetByOverlaps(x = chipPeaks, ranges = up_gene_prom)
dn_ov  <- subsetByOverlaps(x = chipPeaks, ranges = dn_gene_prom)
ctl_ov <- subsetByOverlaps(x = chipPeaks, ranges = ctl_gene_prom)

## Remove shared k27 peaks
shared <- c(subsetByOverlaps(x = up_ov, ranges = dn_ov),
            subsetByOverlaps(x = c(up_ov, dn_ov), ranges = ctl_ov))
up_ov  <- subsetByOverlaps(x = up_ov, ranges = shared, invert = T)
dn_ov  <- subsetByOverlaps(x = dn_ov, ranges = shared, invert = T)
ctl_ov  <- subsetByOverlaps(x = ctl_ov, ranges = shared, invert = T)

## Define function to calculate GC
getGC <- function(x) {
  x %>%
    getSeq(BSgenome.Hsapiens.UCSC.hg38, .) %>%
    letterFrequency(x = ., letters = "GC", as.prob = T) %>%
    as.vector()
}

## Add GC content of promoters as another covariate ##
up_ov$gc <- getGC(up_ov)
dn_ov$gc <- getGC(dn_ov)
ctl_ov$gc <- getGC(ctl_ov)


## Visualize without covariate matching --------------------------------------------------

## Assemble into data.frame
df <- data.frame(
  group = c(
    rep("All peaks in non-diff promoters", length(ctl_ov)),
    rep("Peaks in dn-gene promoters", length(dn_ov)),
    rep("Peaks in up-gene promoters", length(up_ov))
  ),
  peakFC = c(
    ctl_ov$peakFC,
    dn_ov$peakFC,
    up_ov$peakFC
  )
)

## Plot by peakFC
library(ggplot2)
ggplot(data = df, aes(x = group, y = log2(peakFC)))+
  geom_boxplot(outlier.colour = NA) +
  theme_bw()



## Covariate matching & visualization ----------------------------------------------------

## Generate matched control sets
matched <- matchRanges(x = up_ov, univ = ctl_ov, covar = ~peakStrength + gc,
                       method = 'nearest', replace = TRUE)


## Assemble into data.frame
df2 <- data.frame(
  group = c(
    rep("All peaks in non-diff promoters", length(ctl_ov)),
    rep("Matched Controls", length(matched)),
    rep("Peaks in up-gene promoters", length(up_ov))
  ),
  peakFC = c(
    ctl_ov$peakFC,
    matched$peakFC,
    up_ov$peakFC
  )
)

## Plot by peakFC
library(ggplot2)
ggplot(data = df2, aes(x = group, y = log2(peakFC), color = group))+
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(size = 0.01, width = 0.15)+
  theme_bw()

