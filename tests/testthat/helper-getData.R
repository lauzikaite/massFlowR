## use  faahKO package data

files <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)
noise <- 1000
peakwidth <-  c(30, 80)
prefilter <- c(3, 100)
snthresh <- 10
ppm <- 25
integrate <- 1
fitGauss <- FALSE
match <- 1
pearson <- TRUE
