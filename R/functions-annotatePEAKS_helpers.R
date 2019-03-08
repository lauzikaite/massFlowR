annotatePEAKS <- function(dataset = NULL,
                          database = NULL,
                          out_dir = NULL,
                          ncores = 2,
                          rt_err = 2,
                          mz_err = 0.1,
                          n_fts = 10
                          ) {
  if (is.null(dataset)) {
    stop("'dataset' filepath is required")
  }
  if (is.null(database)) {
    stop("'database' filepath is required")
  }
  if (is.null(out_dir)) {
    stop("'out_dir' is required")
  }
  if (!dir.exists(out_dir)) {
    stop("incorrect filepath for 'out_dir' provided")
  }
  ## register paral backend
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  
  ####---- split dataset and database frame into rt regions for parallelisation
  ## order both peak tables by median rt of the peak-groups
  ds_original <- read.csv(dataset, header = T, stringsAsFactors = F)
  db_original <- read.csv(database, header = T, stringsAsFactors = F)
  
  ds_cnames <- c("peakid", "pcs", "mz", "rt", "into")
  ds <- ds_original[, c(match(ds_cnames, colnames(ds_original)))]
  db <- db_original
  
  ## temp fix for float precision
  ds[, c("mz", "rt", "into")] <- t(apply(ds[, c("mz", "rt", "into")], 1, round, digits = 8))
  db[, c("mz", "rt", "into")] <- t(apply(db[, c("mz", "rt", "into")], 1, round, digits = 8))
  
  ds <- orderBYrt(dt = ds, var_name = "pcs")
  db <- orderBYrt(dt = db, var_name = "chemid")
  
  ## get rt/mz error windows
  ds <- addERRS(dt = ds, mz_err = mz_err, rt_err = rt_err)
  db <- addERRS(dt = db, mz_err = mz_err, rt_err = rt_err)
  
  ## get rt region values using ds peak-groups
  ## assign DS peak-groups to bins
  rt_bins <- as.numeric(cut(1:length(unique(ds$pcs)), breaks = ncores))
  for (pcs in unique(ds$pcs)) {
    ds[which(ds$pcs == pcs), "rt_bin"] <- rt_bins[which(unique(ds$pcs) == pcs)]
  }
  ds_bins <- list()
  for (bin in 1:ncores) {
    ds_bins[[bin]] <- ds[which(ds$rt_bin == bin), ]
  }
  
  ## assign DB compounds to bins using DS rt regions
  db_bins <- list()
  for (bin in 1:ncores) {
    rt_val_bin <- min(ds_bins[[bin]]$rt) - rt_err
    rt_val_next <- ifelse(bin < ncores,
                          min(ds_bins[[(bin + 1)]]$rt) - rt_err,
                          Inf)
    db_by_rt <- db[which(db$rt >= rt_val_bin & db$rt < rt_val_next),]
    ## also add peaks that belong to the same chemid
    db_by_cid <- db[which(db$chemid %in% db_by_rt$chemid),]
    db_bins[[bin]] <- db_by_cid
    
  }
  
  ####---- estimate cosines for matching peak-groups between DS and DB
  cos_matches <-
    foreach::foreach(bin = 1:ncores,
                     .inorder = TRUE) %dopar% (
                       massFlowR:::getCOSmat(
                         bin = bin,
                         ds_bin = ds_bins[[bin]],
                         ds_var = "pcs",
                         tmp_bin = db_bins[[bin]],
                         tmp_var = "chemid",
                         mz_err = mz_err,
                         rt_err = rt_err,
                         bins = 0.01
                       )
                     )
  cos_mat <- matrix(0, nrow = length(unique(db$chemid)), ncol = length(unique(ds$pcs)))
  rownames(cos_mat) <- unique(db$chemid)
  colnames(cos_mat) <- unique(ds$pcs)
  
  for(bin in 1:ncores) {
    cos_mat_bin <- cos_matches[[bin]][[1]]
    cos_mat[match(rownames(cos_mat_bin), rownames(cos_mat), nomatch = 0),
            match(colnames(cos_mat_bin), colnames(cos_mat), nomatch = 0)] <-
      cos_mat_bin[match(rownames(cos_mat), rownames(cos_mat_bin), nomatch = 0),
                  match(colnames(cos_mat), colnames(cos_mat_bin), nomatch = 0)]
    
  }
  
  ####---- assign ds peakgroups to db peakgroups using cosines
  cos_assigned <- assignCOS(cos = cos_mat)
  
  ####---- export annotation table
  ds_true <- apply(cos_assigned, 2, function(x) which(x))
  ds_assigned <- which(sapply(ds_true, length) > 0)
  ds_assigned_pcs <- unique(ds$pcs)[ds_assigned]

  db_assigned <- unlist(ds_true[ds_assigned])
  db_assigned_chemid <- unique(db$chemid)[db_assigned]
 
  
  for (x in 1:length(ds_assigned_pcs)) {
    pcs <- ds_assigned_pcs[x]
    chemid <- db_assigned_chemid[x]
    cos <-
      cos_mat[match(chemid, rownames(cos_mat)), match(pcs, colnames(cos_mat))]
    db_chemid <- db[db$chemid == chemid,]
    ds_original[which(ds_original$pcs == pcs), c("chemid",
                                                 "dbid",
                                                 "dbname",
                                                 "cos")] <-
      c(unique(db_chemid[, c("chemid", 
                             "dbid",
                             "dbname"
                             )]),
        cos)
  }
  ## if using full intensity table (temporal fix until decided)
  # ds_cnames <- c(colnames(ds_original)[1:n_fts], c("chemid", "dbid", "dbname", "cos"))
  # ds_cnames_2 <- colnames(ds_original)[!(colnames(ds_original) %in% ds_cnames)]
  # ds_out <- ds_original[ , c(ds_cnames, ds_cnames_2)]
  
  ## if using peak table
  ds_out <- ds_original
  
  write.csv(x = ds_out,
            quote = TRUE,
            # fileEncoding = "utf-8",
            file = file.path(out_dir, "annotated_data.csv"),
            row.names = FALSE)
  
  ## sanity check: do assigned peakgroups pairs have reasonable mz/rt differences
  # for (x in 1:length(ds_assigned_pcs)) {
  # 
  #   pcs <- ds_assigned_pcs[x]
  #   chemid <- db_assigned_chemid[x]
  #   ds_pcs <- ds[ds$pcs == pcs, ]
  #   db_chemid <- db[db$chemid == chemid, ]
  #   if (any(!(min(ds_pcs$mz_l) <= max(db_chemid$mz_h)) |
  #            !(max(ds_pcs$mz_h) >= min(db_chemid$mz_l)))) {
  #     print(x)
  #     break
  #   }
  #   
  #   if (any(!(min(ds_pcs$rt_l) <= max(db_chemid$rt_h)) |
  #           !(max(ds_pcs$rt_h) >= min(db_chemid$rt_l)))) {
  #     print(x)
  #     break
  #   }
  #  # print(cos_mat[match(chemid, rownames(cos_mat)), match(pcs, colnames(cos_mat))])
  # }
  
 
}






orderBYrt <- function(dt, var_name) {
  ## extract unique peak-groups or PCS or chemids
  var_ind <- match(var_name, colnames(dt))
  vars <- unique(dt[, var_ind])
  rt_med <- sapply(vars, function(var) {
    median(dt[which(dt[,var_ind] == var), "rt"])
  })
  var_order <- order(rt_med)
  var_levels <- factor(dt[, var_ind], levels = var_order)
  dt <- dt[order(var_levels),]
  return(dt)
}


getCOSmat <- function(bin, ds_bin, ds_var, tmp_bin, tmp_var, mz_err, rt_err, bins = 0.01) {
  
  ds_var_ind <- match(ds_var, colnames(ds_bin))
  ds_vars <- unique(ds_bin[ , ds_var_ind])
  
  tmp_var_ind <- match(tmp_var, colnames(tmp_bin))
  tmp_vars <- unique(tmp_bin[ , tmp_var_ind])
  
  ## initiate empty matrix to store cosine
  cos_mat <- matrix(0, nrow = length(tmp_vars), ncol = length(ds_vars))
  
  ## initiate empty list to store matching peaks
  mat_list <- vector(mode = "list", length = length(ds_vars))
  
  ## variable is a single DS peak-groups/PCS
  for(var in ds_vars) {
   
    ## get peaks that belong to this variable
    target <- ds_bin[which(ds_bin[, ds_var_ind] == var), ]

    ## match template peaks by mz/rt window
    matches <- apply(target, 1, FUN = getMATCHES, tmp = tmp_bin, tmp_var = tmp_var, target_var = ds_var)
    if (is.null(matches)) {
      next
    }
    matches <- do.call("rbind", matches)
    
    ## build mz spectra using all target peaks and all matched peaks
    matches_all <- tmp_bin[which(tmp_bin[, tmp_var_ind] %in% matches$tmp_var), ]
    all_peaks <- c(target$mz, matches_all$mz)
    mz_spec <-
      data.frame(mz = seq(
        from = min(all_peaks),
        to = max(all_peaks),
        by = bins
      ))
    mz_spec$bin <- 1:nrow(mz_spec)
    
    ## convert target's and tmp peak-groups to vectors
    target_vec <- buildVECTOR(mz_spec = mz_spec, peaks = target)
    
    ## calculate cosines between target peak-group and all matched peak-groups
    ## cos length will be equal to l
    cos <- lapply(unique(matches$tmp_var), function(m) {
      ## extract peaks for the single peak-group
      matched_peaks <- matches_all[which(matches_all[, tmp_var_ind] == m),
                               c("mz", "into")]
      
      ## build vector from matched peaks
      matched_vec <- buildVECTOR(mz_spec = mz_spec, peaks = matched_peaks)

      ## find the cosine of the angle between peakgrs
      cos_angle <-
        (sum(target_vec * matched_vec))  / ((sqrt(sum(
          target_vec * target_vec
        )))  * (sqrt(sum(
          matched_vec * matched_vec
        ))))
      return(cos_angle)
    })
    
    ## update cosine matrix with cosines
    ## rows are tmp variables
    ## columns are ds variables
    cos_mat[match(unique(matches$tmp_var), tmp_vars), which(ds_vars == var)] <- unlist(cos)

    ## update matches list with tmp peaks with which cosine > 0 was obtained
    tmp_vars_pos <- unique(matches$tmp_var)[which(cos > 0)]
    mat_list[[which(ds_vars == var)]] <- matches[which(matches$tmp_var %in% tmp_vars_pos),]
 
  }
  rownames(cos_mat) <- tmp_vars
  colnames(cos_mat) <- ds_vars
  return(list("cos_mat" = cos_mat, "mat_list" = mat_list))
}


getMATCHES <- function(target_peak, tmp, tmp_var, target_var) {
  ## find matching template target_peaks using mz/rt windows of both target_peaks
  mat <- tmp[which((tmp$mz_l >= target_peak["mz_l"] &
                      tmp$mz_l <= target_peak["mz_h"]) |
                     (tmp$mz_h >= target_peak["mz_l"] &
                        tmp$mz_h <= target_peak["mz_h"])),]
  mat <- mat[which((mat$rt_l >= target_peak["rt_l"] &
                      mat$rt_l <= target_peak["rt_h"]) |
                     (mat$rt_h >= target_peak["rt_l"] &
                        mat$rt_h <= target_peak["rt_h"])),]
  if (nrow(mat) > 0) {
    mat$tmp_var <- mat[, match(tmp_var, names(mat))]
    mat$target_var <- target_peak[match(target_var, names(target_peak))]
    mat$target_peakid <- target_peak["peakid"]
    mat <- mat[ ,c("peakid", "mz", "rt", "into", "tmp_var", "target_peakid", "target_var")]
    return(mat)
  } else {
    return(NULL)
  }
}

buildVECTOR <- function(mz_spec, peaks) {
  peaks$bin <- findInterval(x = peaks$mz, mz_spec$mz)
  ## remove duplicating duplicating bins and any bins that are not representing a peak
  peaks <- peaks[which(!duplicated(peaks$bin)), c("bin", "into")]
  mz_spec$into <- 0
  mz_spec[which(mz_spec$bin %in% peaks$bin), "into"] <- peaks$into
  vec <- scaleVEC(mz_spec$into)
  return(vec)
}

assignCOS <- function(cos) {
  
  ## rank by row, i.e. the template vars
  tmp_rank <- t(apply(cos, 1, FUN = rankCOS))
  ## rank by column, i.e. the ds vars
  ds_rank <- apply(cos, 2, FUN = rankCOS)
  ## both matrices have the same dimensions
  
  ## assign peakgroup pairs which have the highest rank both row-wise and column-wise
  assigned <- matrix(FALSE, nrow = nrow(cos), ncol = ncol(cos))
  assigned[which(cos == 0)] <- NA
  assigned[which(ds_rank == tmp_rank & ds_rank == 1)] <- TRUE
  
  while (FALSE %in% assigned) {
    ## add NA to tmp peakgroups that were already assigned with ds peakgroups
    assigned <- apply(assigned, 2, function(x) {
      if (any(na.omit(x))) {
        x[which(x == FALSE)] <- NA
      }
      return(x)
    }
    )
    ## which rows and columns have any un-assigned peakgroups
    tmp_not_assigned <- which(apply(assigned, 1, function(x) {
      if (all(is.na(x))) {
        FALSE
      } else {
        !TRUE %in% x
      }
    }))
    ds_not_assigned <- which(apply(assigned, 2, function(x) {
      if (all(is.na(x))) {
        FALSE
      } else {
        !TRUE %in% x
      }
    }))
    ## for every un-assigned ds peakgroup
    for (var in ds_not_assigned) {
      ## order matches by rank
      current_ds_var <- ds_rank[ , var]
      current_ds_var_ranks <- current_ds_var[order(current_ds_var)]
      for (rank in current_ds_var_ranks) {
        if (rank != 0) {
          current_tmp_var_ind <- which(ds_rank[ , var] == rank)
          ## check if this tmp peak-group is not assigned yet
          ## if this tmp peakgroup is already assigned, assign NA for the pair
          if (current_tmp_var_ind %in% tmp_not_assigned) {
            ## if still free, check if this is the highest rank among remaining ones
            current_tmp_var <- tmp_rank[current_tmp_var_ind, ]
            ## are higher ranks for this tmp var not assigned yet? 
            ## order matches and check every match
            current_tmp_ds_rank <- current_tmp_var[var]
            higher_ranks <- which(current_tmp_var < current_tmp_ds_rank & current_tmp_var != 0)
            higher_ranks <- higher_ranks[order(current_tmp_var[higher_ranks])]
            ## if any of the higher ranking ds are still un-assigned, leave this tmp ds un-assigned now
            if (length(higher_ranks) > 0) {
              if (any(higher_ranks %in% ds_not_assigned)) {
                next
              } else {
                assigned[current_tmp_var_ind, var] <- TRUE
                tmp_not_assigned <- tmp_not_assigned[-c(which(tmp_not_assigned == current_tmp_var_ind))]
                ## finish for this ds peakgroup
                assigned[which(assigned[, var] == FALSE), var] <- NA
                break
              }
            } else {
              ## if this is the highest match for this tmp var
              assigned[current_tmp_var_ind, var] <- TRUE
              tmp_not_assigned <- tmp_not_assigned[-c(which(tmp_not_assigned == current_tmp_var_ind))]
              ## finish for this ds peakgroup
              assigned[which(assigned[, var] == FALSE), var] <- NA
              break
            }
          } else {
            assigned[current_tmp_var_ind, var] <- NA
          }
        }
      }
    }
  }
  return(assigned)
}
  

rankCOS <- function(x) {
  cos_0 <- length(which(x == 0))
  cos_pos <- length(x) - cos_0
  if (cos_pos > 0) {
    cos_ranks <- rep(0, length(x))
    cos_ranks[order(x, decreasing = T)] <- c(
      seq(from = 1, to = cos_pos, by = 1),
      rep(0, cos_0)
    )
  } else {
    cos_ranks <- rep(0, cos_0)
  }
  return(cos_ranks)
}

                   
                    
                    
  
  


