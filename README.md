
# massflowR

[![Build
Status](https://travis-ci.org/lauzikaite/massFlowR.svg?branch=master)](https://travis-ci.org/lauzikaite/massFlowR)

Package for pre-processing of high-resolution, untargeted, centroid
LC-MS data.

`massFlowR` detects and aligns structurally-related spectral peaks
across LC-MS experiment samples.

-----

## Individual samples processing

Each LC-MS file in the experiment is processed independently.

  - Chromatographic peak detection enabled by the *centWave* algorithm
    from `xcms` package.

  - Grouping structurally-related chromatographic peaks into **pseudo
    chemical spectra** (PCS).

Individual samples can be processed in real-time during LC-MS spectra
acquisition, or in parallel, using backend provided by `doParallel`
package.

## Peak alignment

To align chromatographic peaks across all samples in LC-MS experiment,
an algorithm, which compares the overall similarity of pseudo chemical
spectra is implemented (see [Peak
alignment](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/processing.html)).

## Post-alignment processing

Once pseudo chemical spectra are aligned across samples, the obtained
features are validated. Intensity values for each feature in a PCS are
correlated across all samples to identify sets of features with similar
correlation pattern.

Final step in the pipeline is to re-integrate intensity values for
chromatographic peaks that were not detected by the *centWave* using raw
LC-MS files. *m/z* and *rt* integration regions are estimated for each
sample separetely through cubic smoothing spline interpolation.

In the generated feature table, each feature is assigned to a pseudo
chemical spectra, which can be used in manual metabolite annotations.

## Automatic annotation

If in-house chemical reference database is available, obtained PCS are
annotated. For more details how to build a database file, see
[annotation using
database](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/annotation.html)).

-----

## Vignettes

More information is available in vignettes:

  - [massFlowR
    overview](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/massFlowR.html)
  - [Data
    processing](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/processing.html)
  - [Automatic
    annotation](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/vignettes/annotation.html)

-----

# Installation

``` r
# To install dependencies
install_dependencies <- function () {
  dependencies <- c("xcms",  "MSnbase", "faahKO", "igraph", "doParallel",  "foreach", "ggplot2", "viridis", "gridExtra")
  installed <- installed.packages()
  to_install <- subset(dependencies, !(dependencies %in% installed[, "Package"]))
  if (length(to_install) != 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    message("Installing packages: ", to_install, " ...")
    BiocManager::install(to_install)
  } else {
    message("All dependencies are already installed.")
  }
}
install_dependencies()

# To install development version
devtools::install_github("lauzikaite/massflowR", dependencies = FALSE)
```
