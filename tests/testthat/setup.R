test_files <- system.file(c('cdf/WT/wt15.CDF', 'cdf/WT/wt16.CDF'), package = "faahKO")
test_file <- test_files[1]
test_fname <- strsplit(basename(test_files), split = "[.]")[[1]][1]
data_dir <- system.file("testdata", package = "massFlowR")

####---- single datafile preparation
cwt <- xcms::CentWaveParam(ppm = 25,
                           snthresh = 10,
                           noise = 1000,
                           prefilter =  c(3, 100),
                           peakwidth = c(30, 80),
                           integrate = 1,
                           fitgauss = FALSE,
                           verboseColumns = TRUE)
test_raw <-  MSnbase::readMSData(files = test_file, mode = "onDisk")
test_chrom <- xcms::findChromPeaks(object = test_raw, param = cwt)
test_pks <- data.frame(xcms::chromPeaks(test_chrom))
test_pks_rd <- test_pks %>%
  arrange(desc(.data$into)) %>% ## arrange by peak intensity and give a peak number
  mutate(peakid = row_number()) %>%
  group_by(.data$rt, .data$mz) %>%
  arrange(.data$peakid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(peakid = row_number()) %>% ## update peak number after removal of artefactural, duplicating peaks
  data.frame()
test_eic_rd <- xcms::chromatogram(test_raw,
                                    rt = data.frame(
                                      rt_lower = test_pks_rd$rtmin,
                                      rt_upper = test_pks_rd$rtmax),
                                    mz = data.frame(
                                      mz_lower = test_pks_rd$mzmin,
                                      mz_upper = test_pks_rd$mzmax))
test_eic_rd <- lapply(1:nrow(test_eic_rd), function(ch) {
  clean(test_eic_rd[ch, ], na.rm = T)})

####---- multiple datafile prepation with the pipeline
groupPEAKS(files = test_files, out_dir = data_dir, cwt = cwt)

## prepare csv for file input
grouped_files <- list.files(pattern = "peakgrs.csv", path = data_dir, full.names = T)
experiment <- data.frame(filepaths = grouped_files, run_order = 1:2, stringsAsFactors = F)
write.csv(experiment, file.path(data_dir, "experiment.csv"), quote = F, row.names = FALSE)
experiment_file <- file.path(data_dir, "experiment.csv")

## list db template
db_fname <- file.path(data_dir, "DBtemplate.csv")
