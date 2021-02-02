## Load H3K27Ac Peak Data

## Load required libraries
library(data.table)
library(plyranges)
library(magrittr)

## Read in merged peak counts
chipPeaks <-
  fread("inst/extdata/chip/H3K27ac/peaks/h3k27ac_hg38_counts.txt") %>%
  as_granges(keep_mcols = T) %>%
  mutate(start = start + 1)

## Simplify column names
mcols(chipPeaks) %<>%
  set_colnames(gsub('.*WT_(.*)_S_NA_([0-9]).*', '\\1_\\2', colnames(.)))

## Add seqinfo
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqinfo(chipPeaks) <- seqinfo(txdb)[seqnames(seqinfo(chipPeaks))]

## Quantile normalize peak counts
mcols(chipPeaks) %<>%
  as.matrix() %>%
  preprocessCore::normalize.quantiles() %>%
  `colnames<-`(colnames(mcols(chipPeaks)))

## Calculate fold change (PMA/NON) in chip peaks
chipPeaks$peakFC <- 
  rowMedians(as.matrix(mcols(chipPeaks)[3:4])) /
  rowMedians(as.matrix(mcols(chipPeaks)[1:2]))

## Calculate peak strength across samples
chipPeaks$peakStrength <-
  chipPeaks %>%
  plyranges::select(1:4) %>%
  mcols() %>%
  as.matrix %>%
  rowSums()

## Save results to a file
save(chipPeaks, file = "data/h3k27ac_peak_counts.rda")