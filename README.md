
massflowR
=========

Package for pre-processing of high-resolution, untargeted, centroid LC-MS data.

Overview
========

`massflowR` implements pre-processing in two major steps:

1.  Pre-processing of individual datafiles - either in real-time of acquired datafiles one-by-one, or in parallel using backend provided by the `BiocParallel` package:

-   Datafiles are pick-peaked using the *centWave* algorithm \[package `xcms`\].
-   Picked peaks are grouped into **components**.
-   Components are grouped into **clusters**.

1.  Grouping of picked peaks by sample run order

Installation
============

``` r
# For devel version
devtools::install_github("lauzikaite/massflowR")
```
