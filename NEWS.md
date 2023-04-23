# nullranges 1.5.19

* Remove speedglm dependency as it was removed from CRAN
  (April 2023).

# nullranges 1.1.13

* Bug fix that allows nearest matching method to handle
  categorical (or categorical-like) covariates.

# nullranges 1.1.4

* Needed to drop features that have 0 width after trimming in
  bootRanges.
* Made the validity test for bootRanges only look for iter.

# nullranges 1.1.1

* Change to `factor-Rle` output for `bootRanges` to simplify
  the downstream plyranges.

# nullranges 1.0.0

* `nullranges` is released on Bioconductor! the package offers
  the creation of null genomic feature sets, either through
  sampling from a pool in order to match covariates with a 
  particular focal set, or via block bootstrapping of 
  features optionally with respect to a genome segmentation.
  Critically, nullranges is designed as a modular package,
  solely for the purpose of generating null feature sets, 
  and to be used in conjunction with another package for
  calculating overlaps, such as `GenomicRanges` or `plyranges`.
  Let us know your comments, suggestions or feedback on
  Bioconductor support site or through GitHub Issues.
