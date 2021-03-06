% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootRanges.R
\name{bootRanges}
\alias{bootRanges}
\title{Block bootstrap genomic ranges}
\usage{
bootRanges(
  x,
  seg = NULL,
  blockLength,
  R = 1,
  deny = NULL,
  denyOption = c("drop", "trim"),
  proportionLength = TRUE,
  type = c("bootstrap", "permute"),
  withinChrom = FALSE
)
}
\arguments{
\item{x}{the input GRanges}

\item{seg}{the segmentation GRanges, with a column ("state")
indicating segmentation state (optional)}

\item{blockLength}{the length of the blocks (for proportional blocks, this
is the maximal length of block)}

\item{R}{the number of bootstrap samples to generate}

\item{deny}{the GRanges of deny regions (optional)}

\item{denyOption}{whether to "drop" or "trim" bootstrap
ranges that overlap a deny region}

\item{proportionLength}{for segmented block bootstrap,
whether to use scaled block length, (scaling by the proportion
of the segmentation state out of the total genome length)}

\item{type}{the type of null generation (unsegmented only)}

\item{withinChrom}{whether to re-sample (bootstrap) ranges
across chromosomes (default) or only within chromosomes (unsegmented only)}
}
\value{
bootRanges: a GRangesList of length R with the bootstrapped ranges
}
\description{
Performs a block bootstrap (R times) optionally with respect
to a genome segmentation. Returns a bootRanges object, which
is essentially a GRangesList of length R.
}
