% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions-validPEAKS_helpers.R
\name{corPEAKS}
\alias{corPEAKS}
\title{Correlate two peaks' intensities across all samples.}
\usage{
corPEAKS(pair, pkg_ints, min_samples_n)
}
\arguments{
\item{pair}{\code{matrix} with 1 row and 2 columns: "from" and "to", which indicate peaks' peakids.}

\item{pkg_ints}{\code{data.frame} which must contain columns "peakid", "run_order" and "into", 
indicating intensities of every peak in a peak-group in every sample in which they were detected.}

\item{min_samples_n}{\code{numeric} specifying the minimum number of samples in which peak has to be detected in order to be correlated with other peaks.}
}
\value{
Function returns a \code{numeric} indicating Pearson correlation coefficient between the two peaks of interest.
}
\description{
Function takes a pair of peaks and performs Pearson correlation of their "into" values across the samples in which they were detected.
}
