# annotatePEAKS <- function(dataset = NULL,
#                           database = NULL,
#                           out_dir = NULL,
#                           ncores = 2,
#                           rt_err = 2,
#                           mz_err = 0.1,
#                           n_fts = 10
#                           ) {
#   if (is.null(dataset)) {
#     stop("'dataset' filepath is required")
#   }
#   if (is.null(database)) {
#     stop("'database' filepath is required")
#   }
#   if (is.null(out_dir)) {
#     stop("'out_dir' is required")
#   }
#   if (!dir.exists(out_dir)) {
#     stop("incorrect filepath for 'out_dir' provided")
#   }
#   ## register paral backend
#   if (ncores > 1) {
#     cl <- parallel::makeCluster(ncores)
#     doParallel::registerDoParallel(cl)
#   } else {
#     foreach::registerDoSEQ()
#   }
#   
#   ####---- split dataset and database frame into rt regions for parallelisation
#   ## order both peak tables by median rt of the peak-groups
#   ds_original <- read.csv(dataset, header = T, stringsAsFactors = F)
#   db_original <- read.csv(database, header = T, stringsAsFactors = F)
#   
#   ds_cnames <- c("peakid", "pcs", "mz", "rt", "into")
#   ds <- ds_original[, c(match(ds_cnames, colnames(ds_original)))]
#   db <- db_original
#   
#   ## temp fix for float precision
#   ds[, c("mz", "rt", "into")] <- t(apply(ds[, c("mz", "rt", "into")], 1, round, digits = 8))
#   db[, c("mz", "rt", "into")] <- t(apply(db[, c("mz", "rt", "into")], 1, round, digits = 8))
#   
#   ds <- orderBYrt(dt = ds, var_name = "pcs")
#   db <- orderBYrt(dt = db, var_name = "chemid")
#   
#   ## get rt/mz error windows
#   ds <- addERRS(dt = ds, mz_err = mz_err, rt_err = rt_err)
#   db <- addERRS(dt = db, mz_err = mz_err, rt_err = rt_err)
#   
#   ## get rt region values using ds peak-groups
#   ## assign DS peak-groups to bins
#   rt_bins <- as.numeric(cut(1:length(unique(ds$pcs)), breaks = ncores))
#   for (pcs in unique(ds$pcs)) {
#     ds[which(ds$pcs == pcs), "rt_bin"] <- rt_bins[which(unique(ds$pcs) == pcs)]
#   }
#   ds_bins <- list()
#   for (bin in 1:ncores) {
#     ds_bins[[bin]] <- ds[which(ds$rt_bin == bin), ]
#   }
#   
#   ## assign DB compounds to bins using DS rt regions
#   db_bins <- list()
#   for (bin in 1:ncores) {
#     rt_val_bin <- min(ds_bins[[bin]]$rt) - rt_err
#     rt_val_next <- ifelse(bin < ncores,
#                           min(ds_bins[[(bin + 1)]]$rt) - rt_err,
#                           Inf)
#     db_by_rt <- db[which(db$rt >= rt_val_bin & db$rt < rt_val_next),]
#     ## also add peaks that belong to the same chemid
#     db_by_cid <- db[which(db$chemid %in% db_by_rt$chemid),]
#     db_bins[[bin]] <- db_by_cid
#     
#   }
#   
#   ####---- estimate cosines for matching peak-groups between DS and DB
#   cos_matches <-
#     foreach::foreach(bin = 1:ncores,
#                      .inorder = TRUE) %dopar% (
#                        massFlowR:::getCOSmat(
#                          bin = bin,
#                          ds_bin = ds_bins[[bin]],
#                          ds_var = "pcs",
#                          tmp_bin = db_bins[[bin]],
#                          tmp_var = "chemid",
#                          mz_err = mz_err,
#                          rt_err = rt_err,
#                          bins = 0.01
#                        )
#                      )
#   cos_mat <- matrix(0, nrow = length(unique(db$chemid)), ncol = length(unique(ds$pcs)))
#   rownames(cos_mat) <- unique(db$chemid)
#   colnames(cos_mat) <- unique(ds$pcs)
#   
#   for(bin in 1:ncores) {
#     cos_mat_bin <- cos_matches[[bin]][[1]]
#     cos_mat[match(rownames(cos_mat_bin), rownames(cos_mat), nomatch = 0),
#             match(colnames(cos_mat_bin), colnames(cos_mat), nomatch = 0)] <-
#       cos_mat_bin[match(rownames(cos_mat), rownames(cos_mat_bin), nomatch = 0),
#                   match(colnames(cos_mat), colnames(cos_mat_bin), nomatch = 0)]
#     
#   }
#   
#   ####---- assign ds peakgroups to db peakgroups using cosines
#   cos_assigned <- assignCOS(cos = cos_mat)
#   
#   ####---- export annotation table
#   ds_true <- apply(cos_assigned, 2, function(x) which(x))
#   ds_assigned <- which(sapply(ds_true, length) > 0)
#   ds_assigned_pcs <- unique(ds$pcs)[ds_assigned]
# 
#   db_assigned <- unlist(ds_true[ds_assigned])
#   db_assigned_chemid <- unique(db$chemid)[db_assigned]
#  
#   
#   for (x in 1:length(ds_assigned_pcs)) {
#     pcs <- ds_assigned_pcs[x]
#     chemid <- db_assigned_chemid[x]
#     cos <-
#       cos_mat[match(chemid, rownames(cos_mat)), match(pcs, colnames(cos_mat))]
#     db_chemid <- db[db$chemid == chemid,]
#     ds_original[which(ds_original$pcs == pcs), c("chemid",
#                                                  "dbid",
#                                                  "dbname",
#                                                  "cos")] <-
#       c(unique(db_chemid[, c("chemid", 
#                              "dbid",
#                              "dbname"
#                              )]),
#         cos)
#   }
#   ## if using full intensity table (temporal fix until decided)
#   # ds_cnames <- c(colnames(ds_original)[1:n_fts], c("chemid", "dbid", "dbname", "cos"))
#   # ds_cnames_2 <- colnames(ds_original)[!(colnames(ds_original) %in% ds_cnames)]
#   # ds_out <- ds_original[ , c(ds_cnames, ds_cnames_2)]
#   
#   ## if using peak table
#   ds_out <- ds_original
#   
#   write.csv(x = ds_out,
#             quote = TRUE,
#             # fileEncoding = "utf-8",
#             file = file.path(out_dir, "annotated_data.csv"),
#             row.names = FALSE)
#   
#   ## sanity check: do assigned peakgroups pairs have reasonable mz/rt differences
#   # for (x in 1:length(ds_assigned_pcs)) {
#   # 
#   #   pcs <- ds_assigned_pcs[x]
#   #   chemid <- db_assigned_chemid[x]
#   #   ds_pcs <- ds[ds$pcs == pcs, ]
#   #   db_chemid <- db[db$chemid == chemid, ]
#   #   if (any(!(min(ds_pcs$mz_l) <= max(db_chemid$mz_h)) |
#   #            !(max(ds_pcs$mz_h) >= min(db_chemid$mz_l)))) {
#   #     print(x)
#   #     break
#   #   }
#   #   
#   #   if (any(!(min(ds_pcs$rt_l) <= max(db_chemid$rt_h)) |
#   #           !(max(ds_pcs$rt_h) >= min(db_chemid$rt_l)))) {
#   #     print(x)
#   #     break
#   #   }
#   #  # print(cos_mat[match(chemid, rownames(cos_mat)), match(pcs, colnames(cos_mat))])
#   # }
#   
#  
# }
