# checkFILE -------------------------------------------------------------------------------------------------------
## Check and load file for alignment
checkFILE <- function(file = NULL) {
  if (!file.exists(file)) {
    stop("incorrect filepath for: ", file)
  }
  dt <- read.csv(file, stringsAsFactors = F)
  required_colnames <-
    c("peakid",
      "mz",
      "mzmin",
      "mzmax",
      "rt",
      "rtmin",
      "rtmax",
      "into",
      "peakgr")
  if (any(!required_colnames %in% colnames(dt))) {
    stop(
      "incorrect file: ",
      file,
      " \n ",
      "missing columns: ",
      paste0(required_colnames[which(!required_colnames %in% colnames(dt))], collapse = ", ")
    )
  }
  return(dt)
}

# do_alignPEAKS -------------------------------------------------------------------------------------------------------
do_alignPEAKS <- function(ds,
                          tmp,
                          ds_var_name,
                          tmp_var_name,
                          mz_err,
                          rt_err,
                          bins,
                          ncores) {
  
  ####---- split dataset and template into rt regions for parallelisation
  rt_bins <- getRTbins(ds = ds,
                       tmp = tmp,
                       ds_var_name = ds_var_name,
                       tmp_var_name = tmp_var_name,
                       mz_err = mz_err,
                       rt_err = rt_err,
                       ncores = ncores)
  ds_bins <- rt_bins$ds
  tmp_bins <- rt_bins$tmp
  
  ds_vars_by_bins <- lapply(ds_bins, function(x) unique(x$peakgr))
  ds_vars <-  unique(unlist(lapply(ds_bins, "[[", "peakgr")))
  tmp_vars <-  unique(unlist(lapply(tmp_bins, "[[", "peakgr")))
  
  ####---- estimate cosines for matching peak-groups between dataset and template
  if (ncores > 1) {
    cos_matches <-
      foreach::foreach(bin = 1:ncores,
                       .inorder = TRUE) %dopar% (
                         massFlowR:::getCOSmat(
                           bin = bin,
                           ds_bin = ds_bins[[bin]],
                           ds_var = "peakgr",
                           tmp_bin = tmp_bins[[bin]],
                           tmp_var = "peakgr",
                           mz_err = params$mz_err,
                           rt_err = params$rt_err,
                           bins = 0.01
                         )
                       )
  } else {
    ## run with lapply to have the same, 2-list-within-a-bin-list, structure as with foreach
    cos_matches <- lapply(
      X = 1,
      FUN = massFlowR:::getCOSmat,
      ds_bin = ds_bins[[1]],
      ds_var = "peakgr",
      tmp_bin = tmp_bins[[1]],
      tmp_var = "peakgr",
      mz_err = params$mz_err,
      rt_err = params$rt_err,
      bins = 0.01)
  }
  cos_mat <- matrix(0, nrow = length(tmp_vars), ncol = length(ds_vars))
  rownames(cos_mat) <- tmp_vars
  colnames(cos_mat) <- ds_vars
  
  for (bin in 1:ncores) {
    cos_mat_bin <- cos_matches[[bin]][[1]]
    cos_mat[match(rownames(cos_mat_bin), rownames(cos_mat), nomatch = 0),
            match(colnames(cos_mat_bin), colnames(cos_mat), nomatch = 0)] <-
      cos_mat_bin[match(rownames(cos_mat), rownames(cos_mat_bin), nomatch = 0),
                  match(colnames(cos_mat), colnames(cos_mat_bin), nomatch = 0)]
  }
  
  ####---- assign ds peakgroups to tmp according to cosines
  cos_assigned <- assignCOS(cos = cos_mat)
  ds_true <- apply(cos_assigned, 2, function(x) which(x))
  ds_assigned <- which(sapply(ds_true, length) > 0)
  ds_vars_assigned <- ds_vars[ds_assigned]
  tmp_assigned <- unlist(ds_true[ds_assigned])
  tmp_vars_assigned <- tmp_vars[tmp_assigned]
  
  ####---- extract matches for assigned peakgroups only
  ds_to_tmp <- rep(list(NA), length(ds_vars))

  for (var in 1:length(ds_vars)) {
    ds_var <- ds_vars[var]
    if (ds_var %in% ds_vars_assigned) {
      ## get the corresponding assigned tmp peagroup
      tmp_var <- tmp_vars_assigned[which(ds_vars_assigned == ds_var)]
      ## extract mathes between the assigned peakgroups
      bin <- which(sapply(ds_vars_by_bins, function(x) ds_var %in% x))
      mat <- cos_matches[[bin]][[2]][[which(ds_vars_by_bins[[bin]] == ds_var)]]
      mat <- mat[which(mat$tmp_var == tmp_var), ]
      ## extract cosine between the assigned peakgroups
      cos <- cos_mat[match(tmp_var, rownames(cos_mat)), match(ds_var, colnames(cos_mat))]
    } else {
      tmp_var <- NULL
      mat <- NULL
      cos <- NULL
    }
    ds_to_tmp[[var]] <- list("ds" = ds_var, "tmp" = tmp_var, "mat" = mat, "cos" = cos)
  }
  return(ds_to_tmp)
}

# getRTbins -------------------------------------------------------------------------------------------------------
getRTbins <- function(ds, tmp, ds_var_name, tmp_var_name, mz_err, rt_err, ncores) {
  ####---- split dataset-of-interest and template into rt regions (bins) for parallelisation
  ## order both peak tables by median rt of the peak-groups
  ds <- orderBYrt(dt = ds, var_name = ds_var_name)
  tmp <- orderBYrt(dt = tmp, var_name = tmp_var_name)
  
  ## get rt/mz error windows
  ds <- addERRS(dt = ds, mz_err = mz_err, rt_err = rt_err)
  tmp <- addERRS(dt = tmp, mz_err = mz_err, rt_err = rt_err)
  
  ## get rt region values using ds peak-groups
  ## assign DS peak-groups to bins
  ds_var_ind <- match(ds_var_name, colnames(ds))
  tmp_var_ind <- match(tmp_var_name, colnames(tmp))
  
  ## assign bin values to ds peak-groups so that each bin has more or less equal number of peak-groups
  ds_vars <- unique(ds[ , ds_var_ind])
  if (ncores > 1) {
    rt_bins <- as.numeric(cut(1:length(ds_vars), breaks = ncores))
  } else {
    rt_bins <- rep(1, length(ds_vars))
  }
  for (var in ds_vars) {
    ds[which(ds[, ds_var_ind] == var), "rt_bin"] <- rt_bins[which(ds_vars == var)]
  }
  ## split ds frame into corresponding bins and save sub-frames to a list
  ds_bins <- list()
  for (bin in 1:ncores) {
    ds_bins[[bin]] <- ds[which(ds$rt_bin == bin), ]
  }
  
  ## split template frame to bins using ds bins rt values
  tmp_bins <- list()
  for (bin in 1:ncores) {
    ## get all peaks that are below current bin max rt value,
    ## and have not been included in the previous bin
    rt_val_max <- max(ds_bins[[bin]]$rt_h)
    peakgrs_previous <- if (bin == 1) {
      0
    } else {
      tmp_bins[[bin - 1]]$peakgr
    }
    tmp_by_rt <- tmp[which(tmp$rt_h <= rt_val_max &
                             !tmp$peakgr %in% peakgrs_previous), ]
    ## also add peaks that belong to the same peak-group as any of the peaks matched by rt
    tmp_by_group <- tmp[which(tmp[ , tmp_var_ind] %in% tmp_by_rt[ , tmp_var_ind]), ]
    tmp_bins[[bin]] <- tmp_by_group
  }
  return(list(ds = ds_bins, tmp = tmp_bins))
}      

# orderBYrt -------------------------------------------------------------------------------------------------------
orderBYrt <- function(dt, var_name) {
  ## extract unique peak-groups or PCS or chemids
  var_ind <- match(var_name, colnames(dt))
  vars <- unique(dt[, var_ind])
  rt_med <- sapply(vars, function(var) {
    median(dt[which(dt[,var_ind] == var), "rt"])
  })
  var_levels <- factor(dt[, var_ind], levels = vars[order(rt_med)])
  dt_ordered <- dt[order(var_levels),]
  return(dt_ordered)
}

# addERRS ---------------------------------------------------------------------------------------------------------
## add error windows using user-defined mz/rt values
addERRS <- function(dt, mz_err, rt_err) {
  dt[, c("mz_l",
         "mz_h",
         "rt_l",
         "rt_h")] <- c(dt$"mz"-mz_err,
                       dt$"mz"+mz_err,
                       dt$"rt"-rt_err,
                       dt$"rt"+rt_err)
  return(dt)
}

# getCOSmat ---------------------------------------------------------------------------------------------------------
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
    spec <- data.frame(
      mz = seq(
        from = min(all_peaks),
        to = max(all_peaks),
        by = bins
      ))
    spec$bin <- 1:nrow(spec)
    
    ## convert target's and tmp peak-groups to vectors
    target_vec <- buildVECTOR(spec = spec, peaks = target)
    
    ## calculate cosines between target peak-group and all matched peak-groups
    ## cos length will be equal to l
    cos <- lapply(unique(matches$tmp_var), function(m) {
      ## extract peaks for the single peak-group
      matched_peaks <- matches_all[which(matches_all[, tmp_var_ind] == m),
                                   c("mz", "into")]
      
      ## build vector from matched peaks
      matched_vec <- buildVECTOR(spec = spec, peaks = matched_peaks)
      
      ## find the cosine of the angle between peakgrs
      cos_angle <-
        (sum(target_vec * matched_vec))  / ((sqrt(sum(
          target_vec * target_vec
        )))  * (sqrt(sum(
          matched_vec * matched_vec
        ))))
      return(round(cos_angle, 4))
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

# getMATCHES ---------------------------------------------------------------------------------------------------------
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

# buildVECTOR ---------------------------------------------------------------------------------------------------------
buildVECTOR <- function(spec, peaks) {
  peaks$bin <- findInterval(x = peaks$mz, spec$mz)
  ## remove duplicating duplicating bins and any bins that are not representing a peak
  peaks <- peaks[which(!duplicated(peaks$bin)), c("bin", "into")]
  spec$into <- 0
  spec[match(peaks$bin, spec$bin), "into"] <- peaks$into
  vector <- scaleSPEC(spec = spec)
  return(vector)
}

# scaleSPEC ---------------------------------------------------------------------------------------------------------
scaleSPEC <- function(spec, m = 0.6, n = 3) {
  ## Version A - scale to unit length
  spec$into / (sqrt(sum(spec$into * spec$into)))
  
  ## Version B - according to Stein & Scott, 1994
  ## scale the intensity of every mz ion separately
  ## use m and n weighting factors taken from Stein & Scott, 1994
  # apply(spec, 1, function(x) {
  #   x[["into"]] ^ m * x[["mz"]] ^ n
  # })
}

# assignCOS ---------------------------------------------------------------------------------------------------------
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

# rankCOS ---------------------------------------------------------------------------------------------------------
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