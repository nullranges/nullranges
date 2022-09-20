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
while detailed vignettes on matching or bootstrapping can be found
under `Articles`. The `Reference` tab lists function help pages.

## Installation

This package can be installed via Bioconductor:

```
BiocManager::install("nullranges")
```

## Manuscripts

*matchRanges* manuscript:

Eric S. Davis, Wancen Mu, Stuart Lee, Mikhail G. Dozmorov,
Michael I. Love, Douglas H. Phanstiel. (2022)
"matchRanges: Generating null hypothesis genomic ranges
via covariate-matched sampling."
*bioRxiv*
[doi: 10.1101/2022.08.05.502985](https://doi.org/10.1101/2022.08.05.502985)

*bootRanges* manuscript:

Wancen Mu, Eric S. Davis, Stuart Lee, Mikhail G. Dozmorov,
Douglas H. Phanstiel, Michael I. Love.
(2022) "bootRanges: Flexible generation of null sets
of genomic ranges for hypothesis testing."
*bioRxiv*
[doi: 10.1101/2022.09.02.506382](https://doi.org/10.1101/2022.09.02.506382)

## Tidy Ranges Tutorial

Additional tutorial material for performing
[tidy ranges](https://nullranges.github.io/tidy-ranges-tutorial)
analysis is currently being developed.

## Funding

This work was funded by the Chan Zuckerberg Initiative as part of the
EOSS grants.

![](man/figures/czi.png)
