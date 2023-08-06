# nullranges <img id="nullranges_logo" src="man/figures/logo.png" align="right" width="125"/>

[![R build status](https://github.com/nullranges/nullranges/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/nullranges/nullranges/actions/workflows/check-bioc.yml)

## Generation of null ranges via bootstrapping or covariate matching

Modular package for generation of sets of genomic features
representing the null hypothesis. These can take the form of block
bootstrap samples of ranges using the framework of Bickel et al 2010,
or sets of control ranges that are matched across one or more
covariates with a focal set. *nullranges* is designed to be
inter-operable with other packages for analysis of genomic overlap
enrichment, including the *plyranges* Bioconductor package.

An overview vignette can be found at the `Get started` tab above,
including a **decision tree** informing which type of methods may be
most appropriate, whether *matching* or *bootstrapping*.

Detailed vignettes on matching or bootstrapping can be found under
`Articles`. The `Reference` tab lists function help pages. 

## Installation

This package can be installed via Bioconductor:

```
BiocManager::install("nullranges")
```

Installing nullranges and all of its dependencies on Mac or Windows
with binaries will be very fast (a minute or so).

For installing on Ubuntu, note that many *nullranges* packages are
available as binaries, which greatly speeds up installation.
Follow the instructions for 
[r2u](https://eddelbuettel.github.io/r2u/), 
then install the following via `apt`:

```
r-cran-tidyverse r-cran-ks r-cran-speedglm r-cran-data.table 
r-cran-progress r-cran-ggridges r-cran-biocmanager 
r-bioc-rtracklayer r-bioc-genomicalignments r-bioc-interactionset
```

## Papers

*matchRanges* paper:

Eric S. Davis, Wancen Mu, Stuart Lee, Mikhail G. Dozmorov,
Michael I. Love, Douglas H. Phanstiel. (2023)
"matchRanges: Generating null hypothesis genomic ranges
via covariate-matched sampling."
*Bioinformatics*
[doi: 10.1093/bioinformatics/btad197](https://doi.org/10.1093/bioinformatics/btad197)

*bootRanges* paper:

Wancen Mu, Eric S. Davis, Stuart Lee, Mikhail G. Dozmorov,
Douglas H. Phanstiel, Michael I. Love.
(2023) "bootRanges: Flexible generation of null sets
of genomic ranges for hypothesis testing."
*Bioinformatics*
[doi: 10.1093/bioinformatics/btad190](https://doi.org/10.1093/bioinformatics/btad190)

## Tidy Ranges Tutorial

Additional tutorial material for performing
[tidy ranges](https://tidyomics.github.io/tidy-ranges-tutorial)
analysis is currently being developed.

## Funding

This work was funded by the Chan Zuckerberg Initiative as part of the
EOSS grants.

![](man/figures/czi.png)
