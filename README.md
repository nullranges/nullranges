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

You can browse the code or submit an Issue at the following link:

<https://github.com/nullranges/nullranges>

## Installation

This package can be installed via Bioconductor:

```
BiocManager::install("nullranges")
```

## Tidy Ranges Tutorial

Additional tutorial material for performing
[tidy ranges](https://nullranges.github.io/tidy-ranges-tutorial)
analysis is currently being developed.

## Funding

This work was funded by the Chan Zuckerberg Initiative as part of the
EOSS grants.

![CZI](man/figures/czi.png)
