####---- setup -----
dir <- "/Users/el1514/Documents/Projects/massFlowR/packageDEV/"
# files <- read.table("/Volumes/WIN10_Tier2/mzML_100/airwave1/SLNEG/SR_mac.txt", sep = "\t", stringsAsFactors = F)[,1]
files <- read.table("/Volumes/WIN10_Tier2/mzML_100/airwave1/URNEG/SR_mac.txt", sep = "\t", stringsAsFactors = F)[,1]

####---- params ----
## SLNEG
noise <- 600
prefilter <- c(10, 3000)
peakwidth <- c(3, 10)
snthresh <- 5
ppm <- 25
integrate <- 1
verbose <- TRUE
fitGauss <- FALSE
match <- 1
Pearson <- TRUE

## URNEG
noise <- 600
prefilter <- c(10, 10000)
peakwidth <- c(2, 5)
snthresh <- 5
ppm <- 25
integrate <- 1
verbose <- TRUE
fitGauss <- FALSE
match <- 1
Pearson <- TRUE

####----- universal setup ----
f <- files[1]
out_dir <- paste0(dir, basename(dirname(files[1])), "/")
dir.create(out_dir)
fname <- gsub(".mzML", "", basename(f))
out_dir_fname <- ifelse(Pearson == TRUE,  paste0(out_dir, fname, "_Pearson/"), paste0(out_dir, fname, "_Spearman/"))
dir.create(out_dir_fname)

####---- load and peak-pick the sample ----
raw <- MSnbase::readMSData(f, mode = "onDisk")
pks <- pickPEAKS(raw = raw, ppm = ppm, snthresh = snthresh, noise = noise, prefilter = prefilter, peakwidth = peakwidth, integrate = integrate, verbose = verbose, fitGauss = fitGauss)

####---- extract EICs for picked-peaks ----
eic <- extractEIC(raw = raw, pks = pks)

save.image(paste0(out_dir_fname, "pks-eic.RData"))
load("/Users/el1514/Documents/Projects/massFlowR/packageDEV/SLNEG/AIRWAVE_LNEG_ToF06_P58W12_SR_Pearson/pks-eic.RData")
load("/Users/el1514/Documents/Projects/massFlowR/packageDEV/URNEG/AIRWAVE_RNEG_ToF03_U1W12_SR_Pearson/pks-eic.RData")
####---- correlate EIC and build componenents ----
comps <- buildCOMPS(Pearson = Pearson, match = match, thr = 0.95, plot = FALSE, pks = pks, clean = T) # 3.312262 mins
comps <- buildCOMPS(Pearson = Pearson, match = match, thr = 0.95, plot = T, pks = pks, clean = T)

save.image(paste0(out_dir_fname, "pks-eic-comps.RData"))



####---- before loading package, update documentation (icl NAMESPACE)
devtools::document()
devtools::load_all(".")

########################
####----- testibng
usethis::use_testthat()

## coverage (trevis still needs to be setup)
devtools::test()
covr::report()
