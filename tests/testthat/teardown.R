## delete directory with written output
unlink(list.files(pattern = "peakgrs|_aligned|experiment.csv", data_dir, full.names = T), recursive=TRUE)
message("Awesome")
