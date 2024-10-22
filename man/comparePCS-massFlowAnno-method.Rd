% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-massFlowAnno.R
\name{comparePCS,massFlowAnno-method}
\alias{comparePCS,massFlowAnno-method}
\alias{comparePCS}
\title{Compare selected pseudo chemical spectra.}
\usage{
\S4method{comparePCS}{massFlowAnno}(object, pcs = NULL, out_dir = NULL)
}
\arguments{
\item{object}{\code{massFlowAnno} class object.}

\item{pcs}{\code{numeric} specifying two or more pseudo chemical spectra to look at.}

\item{out_dir}{\code{character} specifying desired directory for generated figure. Default set to NULL and figure is plotted to graphical device.}
}
\value{
Method returns three plots comparing selected pseudo chemical spectra.
}
\description{
Method looks into selected pseudo chemical spectra and returns multiple plots:
(1) spectra in the most intense samples; (2) retention time differences; (3) intensity drift during acquisition time.
If out_dir is provided, generated plots are written as png files. Otherwise, plots are plotted on current graphical device.
\code{massFlowAnno} class object must be annotated using \code{\link{annotateDS}} method first.
}
