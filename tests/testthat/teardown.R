## delete directory with written output
unlink(setdiff(
  list.files(data_dir, full.names = T),
  list.files(data_dir, pattern = "DBtemplate.csv", full.names = T)
),
recursive = TRUE)
message("Awesome")
