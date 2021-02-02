## RNA-seq Differential Gene Expression Analysis

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

## Load required libraries
library(data.table)
library(plyranges)
library(tximeta)
library(readr)
library(DESeq2)

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

## Save results to a file
save(de_genes, file = "data/mono_macro_diff_genes.rda")