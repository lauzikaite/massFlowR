####---- sample datafiles
faahko_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
faahko_fname <- strsplit(basename(faahko_file), split = "[.]")[[1]][1]
massFlowR_dir <- file.path(system.file("tests",package="massFlowR"), "objects")
dir.create(massFlowR_dir)

####---- xcms and MSnbase objects
cwt <- xcms::CentWaveParam(ppm = 25,
                                snthresh = 10,
                                noise = 1000,
                                prefilter =  c(3, 100),
                                peakwidth = c(30, 80),
                                integrate = 1,
                                fitgauss = FALSE,
                                verboseColumns = TRUE)
faahko_raw <-  MSnbase::readMSData(files = faahko_file, mode = "onDisk")
faahko_chrom <- xcms::findChromPeaks(object = faahko_raw, param = cwt)
faahko_pks <- data.frame(xcms::chromPeaks(faahko_chrom))
faahko_pks_rd <- faahko_pks %>%
  arrange(desc(into)) %>% ## arrange by peak intensity and give a peak number
  mutate(peakid = row_number()) %>%
  group_by(rt, mz) %>%
  arrange(peakid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(peakid = row_number()) %>% ## update peak number after removal of artefactural, duplicating peaks
  data.frame()
faahko_eic_rd <- xcms::chromatogram(faahko_raw,
                                    rt = data.frame(
                                      rt_lower = faahko_pks_rd$rtmin,
                                      rt_upper = faahko_pks_rd$rtmax),
                                    mz = data.frame(
                                      mz_lower = faahko_pks_rd$mzmin,
                                      mz_upper = faahko_pks_rd$mzmax))
faahko_eic_rd <- lapply(1:nrow(faahko_eic_rd), function(ch) {
  clean(faahko_eic_rd[ch, ], na.rm = T)})

####---- massFlowR objects
data_dir <- system.file("testdata", package = "massFlowR")
studyfiles <- list.files(pattern = "peakgrs.csv", path = data_dir, full.names = T)

## list of 2 different samples
studyfiles <- data.frame(filepaths = studyfiles, run_order = 1:length(studyfiles), stringsAsFactors = F)
write.csv(studyfiles, file.path(data_dir, "studysamples.csv"), quote = F, row.names = FALSE)
study_files <- file.path(data_dir, "studysamples.csv")

## list of 2 identical samples
studyfiles <- data.frame(filepaths = rep(studyfiles$filepaths[1], 2), run_order = 1:2, stringsAsFactors = F)
write.csv(studyfiles, file.path(data_dir, "samestudysamples.csv"), quote = F, row.names = FALSE)
same_files <- file.path(data_dir, "samestudysamples.csv")

## db template
db_fname <- file.path(data_dir, "DBtemplate.csv")


