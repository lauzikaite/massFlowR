## Load objects
faahko_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
faahko_fname <- strsplit(basename(faahko_file), split = "[.]")[[1]][1]
massFlowR_dir <- file.path(system.file(package = "massFlowR"), "tests/objects")
dir.create(massFlowR_dir)

## Run basic object preparation for testing
paramCWT <- xcms::CentWaveParam(ppm = 25,
                                snthresh = 10,
                                noise = 1000,
                                prefilter =  c(3, 100),
                                peakwidth = c(30, 80),
                                integrate = 1,
                                fitgauss = FALSE,
                                verboseColumns = TRUE)
faahko_raw <-  MSnbase::readMSData(files = faahko_file, mode = "onDisk")
faahko_chrom <- xcms::findChromPeaks(object = faahko_raw, param = paramCWT)
faahko_pks <- data.frame(xcms::chromPeaks(faahko_chrom))
faahko_pks_rd <- faahko_pks %>%
  arrange(desc(into)) %>% ## arrange by peak intensity and give a peak number ('pno')
  mutate(peakid = row_number()) %>%
  group_by(rt, mz) %>%
  arrange(peakid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(peakid = row_number()) %>% ## update peak number after removal of duplicating peaks
  data.frame()
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
