## Load libraries ####
library(data.table)
library(tximeta)
library(DESeq2)
library(magrittr)
library(plyranges)

## RNA-seq Differential Gene Expression Analysis -----------------------------------------

## Download and untar Salmon quantification to inst/extdata/rna ##
## The Salmon quant directories are available as a gzipped tar file on Zenodo
## https://zenodo.org/record/4484319#.YBdM2NZOmEI
# 
# ## Download
# download.file(url = "https://zenodo.org/record/4484319/files/mono_macro_diff.tgz?download=1",
#               destfile = "inst/extdata/rna/mono_macro_diff.tgz")
# 
# ## Untar
# untar(tarfile = "inst/extdata/rna/mono_macro_diff.tgz",
#       exdir = "inst/extdata/rna/")


## Read in sample sheet as coldata
coldata <-
  fread("inst/extdata/rna/samplesheet.tsv") %>%
  as.data.frame()

## Check that quant paths are correct
file.exists(coldata$files)

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)

## Add gene symbols
library(org.Hs.eg.db)
gse <- addIds(gse, column = "SYMBOL", gene = T)

## Convert to factors
colData(gse)$Bio_Rep <- as.factor(colData(gse)$Bio_Rep)
colData(gse)$Condition <- as.factor(colData(gse)$Condition)

## Build DESeq object
dds <- DESeqDataSet(gse, design = ~Bio_Rep + Condition)

## Filter out lowly expressed genes (at least 10 counts in at least 4 samples)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]

## Fit model
dds <- DESeq(dds)

## Get DE gene results as GRanges
de_genes <- results(dds,
                    contrast = c("Condition", "PMA", "NON"), 
                    format = "GRanges") %>%
  names_to_column("gene_id") %>%
  filter(padj < 0.01)

## Add gene symbols to de_genes
de_genes$gene_symbol <-
  rowData(gse)$SYMBOL[match(de_genes$gene_id, rowData(gse)$gene_id)]


## Load H3K27Ac Peak Data ----------------------------------------------------------------

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


## Do promoters of upregulated genes also show a directional change in K27ac signal? -----

## Define up/down-regulated genes (will reverse)
up_genes <- de_genes %>% filter(log2FoldChange >= 2)
dn_genes <- de_genes %>% filter(log2FoldChange <= -2)

## Get promoters of up/down-regulated genes
up_gene_prom <- promoters(x = up_genes, upstream = 2000, downstream = 200)
dn_gene_prom <- promoters(x = dn_genes, upstream = 2000, downstream = 200)

## Get chipPeaks that overlap promoters of up/down-regulated genes
ov_up <- subsetByOverlaps(x = chipPeaks, ranges = up_gene_prom)
ov_dn <- subsetByOverlaps(x = chipPeaks, ranges = dn_gene_prom)

## Remove shared chipPeaks
shared <- subsetByOverlaps(x = ov_up, ranges = ov_dn)
ov_up  <- subsetByOverlaps(x = ov_up, ranges = shared, invert = T)
ov_dn  <- subsetByOverlaps(x = ov_dn, ranges = shared, invert = T)

## Collect results in a data frame
df <- data.frame(
  group = c(
    rep("All peaks", length(chipPeaks)),
    rep("Peaks in dn-gene promoters", length(ov_dn)),
    rep("Peaks in up-gene promoters", length(ov_up))
  ),
  peakFC = c(
    chipPeaks$peakFC,
    ov_dn$peakFC,
    ov_up$peakFC
  )
)

## Plot by peakFC
library(ggplot2)
ggplot(data = df, aes(x = group, y = log2(peakFC)))+
  geom_boxplot(outlier.colour = NA) +
  theme_bw()


## Are K27ac peaks in up-gene promoters stronger than all peaks? -------------------------

## Plot log2-transformed Peak strengths for all peaks
chipPeaks %>%
  plyranges::select(peakStrength) %>%
  mcols() %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  plot(main = "K27ac Peak Strength")

## Plot log2-transformed Peak strengths for peaks in up-gene promoters
ov_up %>%
  plyranges::select(peakStrength) %>%
  mcols() %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  lines(col = "blue")

## Add legend
legend('topright',
       legend = c("All peaks", "Peaks in up-gene promoters"),
       text.col = c("black", "blue"),
       cex = 0.80,
       bty = 'n')


## Lets control for peak strength by matching it as a covariate --------------------------

## Find peaks in up-gene promoters and other peaks
ov_up <- subsetByOverlaps(x = chipPeaks, ranges = up_gene_prom)
ov_no <- subsetByOverlaps(x = chipPeaks, ranges = up_gene_prom, invert = T)

## Define as either treatment or control groups
ov_up$group <- 1 # treatment
ov_no$group <- 0 # control

## Subset the control group for faster results
set.seed(123)
s <- sample(x = 1:length(ov_no), size = 20000)
ov_no <- ov_no[s]

## Recombine data & convert to data frame
ov_df <- c(ov_up, ov_no) %>% as.data.frame()

## Create matched control object
library(MatchIt)
matchObj <- matchit(data = ov_df,
                    formula = group ~ peakStrength,
                    method = "nearest")

## Extract the matched data
md <- match.data(matchObj)


## Assess the quality of matches ##

## Summary
summary(matchObj, un = F)

## Propensity score
plot(matchObj, type = "jitter", interactive = F)

## Covariate balance qqplot
plot(matchObj, type = "qq", interactive = FALSE,
     which.xs = c("peakStrength"))

## Visualize distributions of peak strength (matched covariate)
md %>%
  filter(group == 0) %>%
  plyranges::select(peakStrength) %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  plot(main = "K27ac Peak Strength")

md %>%
  filter(group == 1) %>%
  plyranges::select(peakStrength) %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  lines(col = "blue")

legend('topright',
       legend = c("Matched Controls", "Peaks in up-gene promoters"),
       text.col = c("black", "blue"),
       cex = 0.70,
       bty = 'n')


## Lets see if our original observations hold up with the matched data -------------------

## Collect results in a data frame
df2 <- data.frame(
  group = c(
    rep("All peaks", length(chipPeaks)),
    rep("Matched controls", length(md$peakFC[md$group == 0])),
    rep("Peaks in up-gene promoters", length(md$peakFC[md$group == 1]))
  ),
  peakFC = c(
    chipPeaks$peakFC,
    md$peakFC[md$group == 0],
    md$peakFC[md$group == 1]
  )
)

## Plot by peakFC
library(ggplot2)
ggplot(data = df2, aes(x = group, y = log2(peakFC)))+
  geom_boxplot(outlier.colour = NA) +
  theme_bw()


## Lets try this with another covariate (GC-content) -------------------------------------

## Get promoters of upregulated genes
up_gene_prom <- promoters(x = up_genes, upstream = 2000, downstream = 200)

## Find peaks in up-gene promoters and other peaks
ov_up <- subsetByOverlaps(x = chipPeaks, ranges = up_gene_prom)
ov_no <- subsetByOverlaps(x = chipPeaks, ranges = up_gene_prom, invert = T)

## Define as either treatment or control groups
ov_up$group <- 1 # treatment
ov_no$group <- 0 # control

## Subset the control group for faster results
set.seed(123)
s <- sample(x = 1:length(ov_no), size = 20000)
ov_no <- ov_no[s]

## Recombine data
ov_gr <- c(ov_up, ov_no)

## Get sequences using biostrings & BSgenome
library(BSgenome.Hsapiens.UCSC.hg38)
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ov_gr)

## Calculate gc content
ov_gr$gc <-
  letterFrequency(x = seqs, letters = "GC", as.prob = T) %>%
  as.vector()

## Convert to data frame
ov_df <- ov_gr %>% as.data.frame()

## Create matched control object
library(MatchIt)
matchObj <- matchit(data = ov_df,
                    formula = group ~ peakStrength + gc,
                    method = "nearest")

## Extract the matched data
md <- match.data(matchObj)


## Assess the quality of matches ##

## Summary
summary(matchObj, un = F)

## Propensity score
plot(matchObj, type = "jitter", interactive = F)

## Covariate balance qqplot
plot(matchObj, type = "qq", interactive = FALSE,
     which.xs = c("peakStrength", "gc"))

par(mfrow = c(1, 2))
## Visualize distributions of peak strength (matched covariate)
md %>%
  filter(group == 1) %>%
  plyranges::select(peakStrength) %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  plot(main = "K27ac Peak Strength", col = "blue")

md %>%
  filter(group == 0) %>%
  plyranges::select(peakStrength) %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  lines(col = "black")

legend('topright',
       legend = c("Matched Controls", "Peaks in up-gene promoters"),
       text.col = c("black", "blue"),
       cex = 0.70,
       bty = 'n')


## Visualize distributions of gc content (matched covariate)
md %>%
  filter(group == 1) %>%
  plyranges::select(gc) %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  plot(main = "GC content", col = "blue")

md %>%
  filter(group == 0) %>%
  plyranges::select(gc) %>%
  as.matrix() %>%
  log2() %>%
  density() %>%
  lines(col = "black")

legend('topright',
       legend = c("Matched Controls", "Peaks in up-gene promoters"),
       text.col = c("black", "blue"),
       cex = 0.70,
       bty = 'n')
par(mfrow = c(1, 1))

## Lets see if our original observations hold up with the new matched data ---------------

## Collect results in a data frame
df3 <- data.frame(
  group = c(
    rep("All peaks", length(chipPeaks)),
    rep("Matched controls", length(md$peakFC[md$group == 0])),
    rep("Peaks in up-gene promoters", length(md$peakFC[md$group == 1]))
  ),
  peakFC = c(
    chipPeaks$peakFC,
    md$peakFC[md$group == 0],
    md$peakFC[md$group == 1]
  )
)

## Plot by peakFC
library(ggplot2)
ggplot(data = df3, aes(x = group, y = log2(peakFC)))+
  geom_boxplot(outlier.colour = NA) +
  theme_bw()


## Are these differences significant?
library(lmtest)
library(sandwich)

fit <- lm(data = md,
          formula = peakFC ~ peakStrength + gc + group,
          weights = weights)

coeftest(fit, vcov. = vcovCL, cluster = ~subclass)

