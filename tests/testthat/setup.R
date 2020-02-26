test_fnames <- system.file(c('cdf/WT/wt15.CDF', 'cdf/WT/wt16.CDF'), package = "faahKO")
test_fname <- test_fnames[1]
test_basenames <- sapply(test_fnames, function(fname) {
  strsplit(basename(fname), split = "[.]")[[1]][1]
}, USE.NAMES = F)
test_basename <- test_basenames[1]
data_dir <- "~/Documents/Projects/massFlowR/packageDEV/unittest"
# data_dir <- file.path(system.file(package = "massFlowR"), "testdata/")

####---- single datafile preparation
cwt <- xcms::CentWaveParam(ppm = 25,
                           snthresh = 10,
                           noise = 1000,
                           prefilter =  c(3, 100),
                           peakwidth = c(30, 80),
                           verboseColumns = TRUE,
                           mzdiff = 0)
test_raw <-  MSnbase::readMSData(files = test_fname, mode = "onDisk")
test_res <- xcms::findChromPeaks(object = test_raw, param = cwt)
test_pks <- data.frame(xcms::chromPeaks(test_res))

## arrange by peak intensity and give a peak number
test_pks_rd <- test_pks[order(test_pks$into, decreasing = T), ]
test_pks_rd$peakid <- 1:nrow(test_pks_rd)
## remove artefactural, duplicating peaks
pks_unique <- unique(test_pks_rd[,c("mz", "rt")])
test_pks_rd <- lapply(1:nrow(pks_unique),
                    FUN = cleanPEAKS,
                    dt_unique = pks_unique,
                    dt = test_pks_rd)
test_pks_rd <- do.call("rbindCLEAN", test_pks_rd)
test_pks_rd$peakid <- 1:nrow(test_pks_rd)


test_eic_rd <- xcms::chromatogram(
  test_raw,
  rt = data.frame(rt_lower = test_pks_rd$rtmin,
                  rt_upper = test_pks_rd$rtmax),
  mz = data.frame(mz_lower = test_pks_rd$mzmin,
                  mz_upper = test_pks_rd$mzmax)
)
test_eic_rd <- lapply(1:nrow(test_eic_rd), function(ch) {
  MSnbase::clean(test_eic_rd[ch, ], na.rm = T)})

## prepare metadata
metadata <- data.frame(filename = test_basenames,
                       is_sr = c(TRUE, TRUE),
                       run_order = 1:2,
                       raw_filepath = test_fnames,
                       proc_filepath = paste0(file.path(data_dir, test_basenames), "_peakgrs.csv"),
                       stringsAsFactors = F
                       )
write.csv(metadata, file.path(data_dir, "metadata.csv"), quote = F, row.names = FALSE)
meta_fname <- file.path(data_dir, "metadata.csv")

####---- multiple datafile prepation with the pipeline
groupPEAKS(file = meta_fname, out_dir = data_dir, cwt = cwt)
grouped_fnames <- paste0(file.path(data_dir, test_basenames), "_peakgrs.csv")

# Duplicated tables -----------------------------------------------------------------------------------------------
## write two duplicated csv for sample wt15
single_table <- read.csv(grouped_fnames[1], header = T, stringsAsFactors = F)
dup_basenames <- c("test_file1", "test_file2")
dup_fnames <- paste0(file.path(data_dir, dup_basenames), "_peakgrs.csv")
write.csv(single_table, dup_fnames[1], quote = F, row.names = FALSE)
write.csv(single_table, dup_fnames[2], quote = F, row.names = FALSE)
metadata_dup <- data.frame(filename = dup_basenames,
                           is_sr = c(TRUE, TRUE),
                           run_order = 1:2, 
                           raw_filepath = rep(test_fnames[1], 2),
                           proc_filepath = dup_fnames
                           )
write.csv(metadata_dup, file.path(data_dir, "metadata_dup.csv"), quote = F, row.names = FALSE)
dup_meta_fname <- file.path(data_dir, "metadata_dup.csv")

# Noisy tables -----------------------------------------------------------------------------------------------------
## write a sligthly modified peak table for sample wt15
## take the biggest peakgroup and add it again with modified peaks
biggest_pkg <- order(table(single_table$peakgr), decreasing = T)[1]

## take only half of the peaks in the biggest peakgr
biggest <- single_table[which(single_table$peakgr == biggest_pkg), ]
biggest <- biggest[1:(nrow(biggest)/2), ]

## make noisy table, which contains all original peakgrs, but the biggest
noisy <- rbind(single_table[single_table$peakgr != biggest_pkg, ],
               ## add (half the peaks - 1) of the biggests peakgr 
               biggest[-nrow(biggest), ])
noisy$noisy[which(noisy$peakgr == biggest_pkg)] <- 1

## modify mz/rt/into values of the biggest peakgr and give new peakgr id
noisy_pkg <- max(single_table$peakgr) + 1

rep_biggest <- biggest
rep_biggest$peakid <- NA # will assign new once added to the table
rep_biggest$into <- rep_biggest$into * 1.5 * 1:nrow(rep_biggest)
rt_new <- min(rep_biggest$rt) - 0.1
rep_biggest$rtmin <- rep_biggest$rtmin - (rep_biggest$rt - rt_new)
rep_biggest$rtmax <- rep_biggest$rtmax - (rep_biggest$rt - rt_new)
rep_biggest$rt <- rt_new
rep_biggest$peakgr <- noisy_pkg
rep_biggest$noisy <- 2

## add modified biggest peakgr as a new peakgr
noisy <- rbind(noisy,
               rep_biggest)

## arrange by peak intensity and give a peak number
noisy <- noisy[order(noisy$into, decreasing = T), ]
noisy$peakid <- 1:nrow(noisy)

noisy_basenames <- c("test_file1", "test_file3")
noisy_fnames <- paste0(file.path(data_dir, noisy_basenames), "_peakgrs.csv")
write.csv(noisy, noisy_fnames[2], quote = F, row.names = FALSE)
noisy_metadata <-
  data.frame(
    filename = noisy_basenames,
    is_sr = c(TRUE, TRUE),
    run_order = 1:2,
    raw_filepath = rep(test_fnames[1], 2),
    proc_filepath = noisy_fnames
  )
write.csv(noisy_metadata, file.path(data_dir, "metadata_noisy.csv"), quote = F, row.names = FALSE)
noisy_meta_fname <- file.path(data_dir, "metadata_noisy.csv")


# large study tables -----------------------------------------------------------------------------------------------------
large_fnames <-
  list.files(file.path(system.file(package = "faahKO"), "cdf/WT/"), full.names = T)
large_basenames <- sapply(large_fnames, function(fname) {
  strsplit(basename(fname), split = "[.]")[[1]][1]
}, USE.NAMES = F)
## prepare metadata
large_metadata <- data.frame(filename = large_basenames,
                             is_sr = c(TRUE, TRUE),
                             run_order = 1:length(large_basenames),
                             raw_filepath = large_fnames,
                             proc_filepath = paste0(file.path(data_dir, large_basenames), "_peakgrs.csv"),
                             stringsAsFactors = F)
write.csv(large_metadata, file.path(data_dir, "metadata_large.csv"), quote = F, row.names = FALSE)
large_meta_fname <- file.path(data_dir, "metadata_large.csv")
groupPEAKS(file = large_meta_fname, out_dir = data_dir, cwt = cwt)

# Real-time implementation ------------------------------------------------


# other vars -----------------------------------------------------------------------------------------------------
## list db template
## template includes three first peak-groups of sample wt15
db_fname <- file.path(data_dir, "DBtemplate.csv")

## alignPEAKS parameters
rt_err <- 10
mz_err <- 0.01
bins <- 0.01
cutoff <- 0.5
min_samples_prop <- 0.3

## prep for the default 2-core implementation for low-level un-exported functions
ncores <- 2
doParallel::registerDoParallel(cores = ncores)
