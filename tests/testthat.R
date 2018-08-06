library(testthat)
library(massflowR)
library(faahKO)

####---- Create objects to be used with multiple tests
## Use faahKO package data
# faahko_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#                   system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#                   system.file('cdf/KO/ko18.CDF', package = "faahKO"))

faahko_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
faahko_fname <- strsplit(basename(faahko_file), split = "[.]")[[1]][1]

## Output directory for txt output
massflowR_dir <- system.file(package = "massflowR")

## centWave parameters and groupCOMPS() parameters for faahKO files
paramCWT <- xcms::CentWaveParam(ppm = 25,
                                snthresh = 10,
                                noise = 1000,
                                prefilter =  c(3, 100),
                                peakwidth = c(30, 80),
                                integrate = 1,
                                fitgauss = FALSE,
                                verboseColumns = TRUE)

####---- XCMS-based objects

## An OnDiskMSnExp object using single file
faahko_raw <-  MSnbase::readMSData(files = faahko_file, mode = "onDisk")

## Peak table
faahko_chrom <- xcms::findChromPeaks(object = faahko_raw, param = paramCWT)
faahko_pks <- data.frame(xcms::chromPeaks(faahko_chrom))
faahko_pks_rd <- faahko_pks %>%
  ## arrange by peak intensity and give a peak number ('pno')
  arrange(desc(into)) %>%
  mutate(pno = row_number()) %>%
  group_by(rt, mz) %>%
  arrange(pno) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  ## update peak number after removal of duplicating peaks
  mutate(pno = row_number()) %>%
  data.frame()

## List of EICs for picked peaks
faahko_eic <- xcms::chromatogram(faahko_raw,
                                 rt = data.frame(
                                   rt_lower = faahko_pks$rtmin,
                                   rt_upper = faahko_pks$rtmax),
                                 mz = data.frame(
                                   mz_lower = faahko_pks$mzmin,
                                   mz_upper = faahko_pks$mzmax))

faahko_eic_rd <- xcms::chromatogram(faahko_raw,
                                    rt = data.frame(
                                      rt_lower = faahko_pks_rd$rtmin,
                                      rt_upper = faahko_pks_rd$rtmax),
                                    mz = data.frame(
                                      mz_lower = faahko_pks_rd$mzmin,
                                      mz_upper = faahko_pks_rd$mzmax))

faahko_eic_rd <- lapply(1:nrow(faahko_eic_rd), function(ch) {
  clean(faahko_eic_rd[ch, ], na.rm = T)})

####---- Initiate testing functions now
test_check("massflowR")


