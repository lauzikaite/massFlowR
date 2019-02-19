test_fnames <- system.file(c('cdf/WT/wt15.CDF', 'cdf/WT/wt16.CDF'), package = "faahKO")
test_fname <- test_fnames[1]
test_basenames <- sapply(test_fnames, function(fname) {
  strsplit(basename(fname), split = "[.]")[[1]][1]
}, USE.NAMES = F)
test_basename <- test_basenames[1]
data_dir <- "~/Documents/Projects/massFlowR/packageDEV/unittest/"

####---- single datafile preparation
cwt <- xcms::CentWaveParam(ppm = 25,
                           snthresh = 10,
                           noise = 1000,
                           prefilter =  c(3, 100),
                           peakwidth = c(30, 80),
                           verboseColumns = TRUE)
test_raw <-  MSnbase::readMSData(files = test_fname, mode = "onDisk")
test_res <- xcms::findChromPeaks(object = test_raw, param = cwt)
test_pks <- data.frame(xcms::chromPeaks(test_res))
test_pks_rd <- test_pks %>%
  ## arrange by peak intensity and give a peak number
  arrange(desc(into)) %>% 
  mutate(peakid = row_number()) %>%
  group_by(rt, mz) %>%
  arrange(peakid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  ## update peak number after removal of artefactural, duplicating peaks
  mutate(peakid = row_number()) %>% 
  data.frame()
test_eic_rd <- xcms::chromatogram(test_raw,
                                    rt = data.frame(
                                      rt_lower = test_pks_rd$rtmin,
                                      rt_upper = test_pks_rd$rtmax),
                                    mz = data.frame(
                                      mz_lower = test_pks_rd$mzmin,
                                      mz_upper = test_pks_rd$mzmax))
test_eic_rd <- lapply(1:nrow(test_eic_rd), function(ch) {
  MSnbase::clean(test_eic_rd[ch, ], na.rm = T)})

####---- multiple datafile prepation with the pipeline
groupPEAKS(files = test_fnames, out_dir = data_dir, cwt = cwt)

## prepare csv for file input
grouped_fnames <- list.files(pattern = "peakgrs.csv", path = data_dir, full.names = T)
metadata <- data.frame(filename = test_basenames,filepaths = grouped_fnames, run_order = 1:2, stringsAsFactors = F)
write.csv(metadata, file.path(data_dir, "metadata.csv"), quote = F, row.names = FALSE)
meta_fname <- file.path(data_dir, "metadata.csv")


# Duplicated tables -----------------------------------------------------------------------------------------------
## write two duplicated csv for sample wt15
single_table <- read.csv(grouped_fnames[1], header = T, stringsAsFactors = F)
dup_fnames <- c(file.path(data_dir, "test_file1.csv"), file.path(data_dir, "test_file2.csv"))
write.csv(single_table, dup_fnames[1], quote = F, row.names = FALSE)
write.csv(single_table, dup_fnames[2], quote = F, row.names = FALSE)
metadata_dup <- data.frame(filename = c("test_file1", "test_file2"), filepaths = dup_fnames, run_order = 1:2, stringsAsFactors = F)
write.csv(metadata_dup, file.path(data_dir, "metadata_dup.csv"), quote = F, row.names = FALSE)
dup_meta_fname <- file.path(data_dir, "metadata_dup.csv")


# Noisy tables -----------------------------------------------------------------------------------------------------
## write a sligthly modified peak table for sample wt15
## take the biggest peakgroup and add it again with modified peaks
biggest_pkg <- order(table(single_table$peakgr), decreasing = T)[1]

## take only half of entries for the biggest peakgr
biggest <- single_table %>%
  filter(peakgr == biggest_pkg) %>%
  slice(1:(n()/2))

noisy <- bind_rows(single_table %>%
                    filter(peakgr != biggest_pkg),
                  ## return the original peakgr with half and minus 1 peaks
                  biggest %>% 
                    slice(1:(n()- 1))) %>% 
  ## modify mz/into values of the replicated peakgr
  bind_rows(biggest %>% 
              ## increase into values so that added peakgr would be higherd in peak table
              mutate(into = into * 1.5 * row_number(),
                     ## modify rt so that added peakgr would be in a different cluster
                     rt_new = min(biggest$rt) - 0.1) %>% 
              mutate(rtmin = rtmin - (rt - rt_new), rtmax = rtmax - (rt - rt_new),
                     rt = rt_new,
                     peakgr = max(single_table$peakgr) + 1,
                     peakid = NULL) %>% 
              select(- rt_new))
noisy_pkg <- max(noisy$peakgr)  
noisy <- noisy %>% 
  ## arrange by peak intensity and give a peak number
  arrange(desc(into)) %>% 
  mutate(peakid = row_number())

noisy_fnames <- c(file.path(data_dir, "test_file1.csv"), file.path(data_dir, "test_file3.csv"))
write.csv(noisy, noisy_fnames[2], quote = F, row.names = FALSE)
noisy_metadata <-
  data.frame(
    filename = c("test_file1", "test_file3"),
    filepaths = noisy_fnames,
    run_order = 1:2,
    stringsAsFactors = F
  )
write.csv(noisy_metadata, file.path(data_dir, "metadata_noisy.csv"), quote = F, row.names = FALSE)
noisy_meta_fname <- file.path(data_dir, "metadata_noisy.csv")

## list db template
## template includes three first peak-groups of sample wt15
db_fname <- file.path(data_dir, "DBtemplate.csv")

## alignPEAKS parameters
rt_err <- 10
mz_err <- 0.01
bins <- 0.01
