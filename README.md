
massflowR
=========

Package for pre-processing of high-resolution, untargeted, centroid LC-MS data.

`massFlowR` detects and aligns structurally-related spectral peaks across LC-MS experiment samples.

------------------------------------------------------------------------

Individual samples processing
-----------------------------

Each LC-MS file in the experiment is processed independently.

-   Chromatographic peak detection enabled by the *centWave* algorithm from `xcms` package.

-   Grouping of chromatographic peaks into structurally-related groups of peaks (see [Peak Grouping](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/massFlowR.html)).

Individual samples can be processed in real-time during LC-MS spectra acquisition, or in parallel, using backend provided by `doParallel` package.

Peak alignment
--------------

To align peaks across all samples in LC-MS experiment, an algorithm, which compares the overall similarity of structurally-related groups of peaks is implemented (see [Peak alignment](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/massFlowR.html)).

Post-alignment processing
-------------------------

Once samples are aligned, the obtained peak groups are validated. Intensity values for each peak in a group are correlated across all samples to identify sets of peaks exhibiting similar behaviour. These sets of peaks are called **Pseudo Chemical Spectra** (PCS).

Final step in the pipeline is to re-integrate intensity for peaks that were not detected by the *centWave* using raw LC-MS files. *m/z* and *rt* values for intensity integration are estimated for each sample separetely. *m/z* and *rt* values are intrapolated using a cubic smoothing spline.

Peak annotation
---------------

If in-house chemical reference database is available, PCS are annotated. For more details how to build a database file, see [annotation using database](https://htmlpreview.github.io/?https://github.com/lauzikaite/massFlowR/blob/master/doc/massFlowR.html)).

------------------------------------------------------------------------

Installation
============

``` r
# For devel version
devtools::install_github("lauzikaite/massflowR")
```
