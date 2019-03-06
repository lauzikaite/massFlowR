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
  mat <- mat[which((mat$mz_l >= target_peak["mz_l"] &
                      mat$mz_l <= target_peak["mz_h"]) |
                     (mat$mz_h >= target_peak["mz_l"] &
                        mat$mz_h <= target_peak["mz_h"])),]
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

                   
                    
                    
  
  


