# show ------------------------------------------------------------------------------------------------------
#' @include classes.R
#'
#' @rdname massFlowAnno-class
#'
#' @param object \code{massFlowAnno} class object
#'
#' @export
#'
setMethod("show", signature = "massFlowAnno", function(object) {
  cat(
    "A \"massFlowAnno\" object with",
    nrow(object@samples),
    " samples and ",
    length(unique(object@db$chemid)),
    "database compounds"
  )
})

# setValidity -----------------------------------------------------------------------------------------------------
setValidity("massFlowAnno", function(object)
  validmassFlowAnno(object))

# annotateDS --------------------------------------------------------------
#' @aliases annotateDS
#'
#' @title Annotate dataset using chemical reference database.
#'
#' @description Method annotates a dataset using a chemical reference database.
#' Peak-groups, or Pseudo Chemical Spectra (PCS), in the dataset are compared with the peak-groups representing reference compounds in the database.
#' PCS are annotated using dot-product estimation between the spectra of the PCS and the corresponding chemical in the database.
#'
#' @param object  \code{massFlowAnno} class object.
#' @param db_file \code{character} for absolute path to the csv file of the chemical reference database, built using \code{\link{buildDB}} function.
#' @param out_dir \code{character} specifying desired directory for output.
#' @param mz_err \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.01.
#' @param rt_err \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 10 (sec).
#' @param ncores \code{numeric} for number of parallel workers to be used. Set 1 for serial implementation. Default set to 2.
#'
#' @return Method writes an updated dataset table with columns 'chemid', 'dbname', and 'similarity',
#' indicating matched chemical reference standards from the database.
#' Column 'similarity' indicates the stringth of spectral similarity between the PCS and the corresponding chemical compound, ranging from 0 to 1.
#'
#' @export
#'
#' @seealso \code{\link{buildDB}}, \code{\link{buildANNO}}
#'
setMethod("annotateDS",
  signature = "massFlowAnno",
  function(object,
             db_file = NULL,
             out_dir = NULL,
             mz_err = 0.01,
             rt_err = 10,
             ncores = 2) {
    ## object validation
    if (!validObject(object)) {
      stop(validObject(object))
    }
    ##  check all inputs
    if (is.null(db_file)) {
      stop("'db_file' is required")
    }
    if (!file.exists(db_file)) {
      stop("incorrect filepath for 'db_file' provided")
    }
    if (is.null(out_dir)) {
      stop("'out_dir' is required!")
    }
    if (!dir.exists(out_dir)) {
      stop("incorrect filepath for 'out_dir' provided")
    }
    if (ncores < 1 | !is.numeric(ncores)) {
      warning("'ncores' was not correctly set. Switching to ncores = 1 (serial performance)")
    }
    ## register paral backend
    if (ncores > 1) {
      doParallel::registerDoParallel(cores = ncores)
    } else {
      foreach::registerDoSEQ()
    }

    ## load and check database table
    db <- read.csv(db_file, header = TRUE, stringsAsFactors = FALSE)
    req_cnames <- c("peakid", "mz", "rt", "into", "chemid", "dbid", "dbname")
    if (any(!req_cnames %in% names(db))) {
      stop("database table must contain columns: ", paste0(req_cnames, collapse = ", "))
    }

    #### ---- compare dataset with the database
    message("Annotating dataset... ")
    ds <- object@ds
    ds_to_db_anno <- do_alignPEAKS(
      ds = ds,
      tmp = db,
      ds_var_name = "pcs",
      tmp_var_name = "chemid",
      mz_err = mz_err,
      rt_err = rt_err,
      bins = 0.1,
      ncores = ncores,
      cutoff = 0,
      anno = TRUE
    )
    feats_mat <- ds_to_db_anno[[2]]
    cos_mat <- ds_to_db_anno[[1]]

    #### ---- prep anotation table for export
    ds_to_db <- lapply(unique(ds$pcs), function(pcs) {
      pcs_mat <- cos_mat[, match(pcs, colnames(cos_mat))]
      pcs_mat <- pcs_mat[order(pcs_mat, decreasing = TRUE)]
      cos <- pcs_mat[which(pcs_mat > 0)]
      chemids <- names(which(pcs_mat > 0))
      dbids <- db$dbid[match(chemids, unique(db$chemid))]
      dbnames <- db$dbname[match(chemids, unique(db$chemid))]
      data.frame(
        pcs = pcs,
        chemid = paste0(chemids, collapse = ", "),
        dbid = paste0(dbids, collapse = ", "),
        dbname = paste0(dbnames, collapse = ", "),
        similarity = paste0(cos, collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
    ds_to_db <- do.call("rbind", ds_to_db)
    write.csv(ds_to_db, file = file.path(out_dir, "DSannotated.csv"), row.names = FALSE)

    object@db <- db
    object@mat <- cos_mat
    object@anno <- feats_mat
    object@params <- list(
      mz_err = mz_err,
      rt_err = rt_err
    )

    if (ncores > 1) {
      foreach::registerDoSEQ()
    }
    if (validObject(object)) {
      message("Dataset was annotated succesfully.")
      return(object)
    }
  }
)


# findANNOchemid ----------------------------------------------------------
#' @aliases findANNOchemid
#' 
#' @title Find annotations for selected chemical compound id
#' 
#' @description Method returns all valid annotations for the selected chemical compound.
#'
#' @param object \code{massFlowAnno} class object.
#' @param chemid \code{numeric} specifying one chemical compound id to look at.
#' @param cutoff \code{numeric} specifying cosine threshold, annotations below are not reported. Default set to 0.
#' 
#' @return Method returns a \code{data.frame} with all valid pseudo chemical spectra annotated to the selected chemical reference database compound.
#' 
#' @export
#'
setMethod("findANNOchemid",
  signature = "massFlowAnno",
  function(object,
             chemid = NULL,
             cutoff = 0) {
    if (is.null(chemid)) {
      stop("'chemid' is required")
    }
    if (length(chemid) > 1) {
      stop("only one 'chemid' is permitted")
    }
    cos_mat <- object@mat
    db <- object@db

    ## extract PCS matching to specified CHEMID
    chemid_mat <- cos_mat[match(chemid, rownames(cos_mat)), ]
    chemid_pcs <- names(which(chemid_mat > cutoff))

    if (length(chemid_pcs) == 0) {
      warning("chemical compound matches 0 PCS with cos > cutoff")
      return(NULL)
    } else {
      ## order matched pcs by their cosines
      cos <- cos_mat[as.character(chemid), chemid_pcs]
      chemid_pcs <- chemid_pcs[order(cos, decreasing = TRUE)]
      cos <- round(cos[order(cos, decreasing = TRUE)], digits = 3)
      res <- data.frame(
        chemid = chemid,
        dbid = unique(db$dbid[db$chemid == chemid]),
        dbname = unique(db$dbname[db$chemid == chemid]),
        pcs = chemid_pcs,
        cos = cos,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
      message(
        "chemical compound matches ", length(chemid_pcs), " PCS with cos > cutoff: \n",
        paste0(chemid_pcs, collapse = ", ")
      )
      return(res)
    }
  }
)

# findANNOpcs -------------------------------------------------------------
#' @aliases findANNOpcs
#' 
#' @title Find annotations for selected pseudo chemical spectra
#' 
#' @description Method returns all valid annotations for the selected pseudo chemical spectra.
#'
#' @param object \code{massFlowAnno} class object.
#' @param pcs \code{numeric} specifying one pseudo chemical spectra to look at.
#' @param cutoff \code{numeric} specifying cosine threshold, annotations below are not reported. Default set to 0.
#'
#' @return Method returns a \code{data.frame} with all valid chemical reference database compounds annotated to the selected pseudo chemical spectra. 
#' 
#' @export
#'
setMethod("findANNOpcs",
  signature = "massFlowAnno",
  function(object,
             pcs = NULL,
             cutoff = 0) {
    if (is.null(pcs)) {
      stop("'pcs' is required")
    }
    if (length(pcs) > 1) {
      stop("only one 'pcs' is permitted")
    }
    cos_mat <- object@mat
    db <- object@db

    ## extract PCS matching to specified CHEMID
    pcs_mat <- cos_mat[, match(pcs, colnames(cos_mat))]
    pcs_chemid <- names(which(pcs_mat > cutoff))

    if (length(pcs_chemid) == 0) {
      warning("PCS matches 0 chemical compounds with cos > cutoff")
      return(NULL)
    } else {
      ## order matched pcs by their cosines
      cos <- cos_mat[pcs_chemid, as.character(pcs)]
      pcs_chemid <- pcs_chemid[order(cos, decreasing = TRUE)]
      cos <- round(cos[order(cos, decreasing = TRUE)], digits = 3)
      res <- data.frame(
        pcs = pcs,
        chemid = pcs_chemid,
        dbid = db$dbid[match(pcs_chemid, db$chemid)],
        dbname = db$dbname[match(pcs_chemid, db$chemid)],
        cos = cos,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
      message(
        "PCS matches ", length(pcs_chemid), " chemical compounds with cos > cutoff: \n",
        paste0(pcs_chemid, collapse = ", ")
      )
      return(res)
    }
  }
)

# checkANNOTATION ---------------------------------------------------------
#' @aliases checkANNOTATION
#'
#' @title Check annotation results between selected database chemical compound and dataset pseudo chemical spectra.
#'
#' @description Method returns summary plots of annotation between the selected chemical compound and pseudo chemical spectra in the dataset.
#' If out_dir is provided, generated plots are written as png files. Otherwise, plots are plotted on current graphical device.
#' \code{massFlowAnno} class object must be annotated using \code{\link{annotateDS}} method first.
#'
#' @param object \code{massFlowAnno} class object.
#' @param chemid \code{numeric} specifying one chemical compound id to look at.
#' @param pcs \code{numeric} specifying one pseudo chemical spectra to look at.
#' @param out_dir \code{character} specifying desired directory for output.
#'
#' @return Method returns two plots and write a csv file: (1) spectra of the selected chemical compound and pseudo chemical spectra from the dataset;
#' (2) retention time of the selected chemical compound and pseudo chemical spectra from the dataset; (3) table with the selected compound and pseudo chemical spectra.
#'
#' @export
#'
setMethod("checkANNOTATION",
  signature = "massFlowAnno",
  function(object,
             chemid = NULL,
             pcs = NULL,
             out_dir = NULL) {
    if (is.null(chemid)) {
      stop("'chemid' is required")
    }
    if (length(chemid) > 1) {
      stop("only one 'chemid' is permitted")
    }
    if (is.null(pcs)) {
      stop("'pcs' is required")
    }
    if (length(pcs) > 1) {
      stop("only one 'pcs' is permitted")
    }
    ## extract features for PCS and CHEMID
    ds <- object@ds
    data_mat <- object@data
    samples_ind <- match(object@samples$filename, colnames(data_mat))
    data_mat <- data_mat[, samples_ind]
    db <- object@db
  
    ## extract features for PCS:
    ## (1) extract intensity values in the most abundant sample
    ## (2) correlate intensities with the main adduct
    pcs_ds <- lapply(pcs, ds = ds, data_mat = data_mat, FUN = prepPCS)
    pcs_ds <- do.call(rbind, pcs_ds)
    
    ## extract features from CHEMID
    chemid_db <- db[db$chemid == chemid, ]
    chemid_db$into_scaled <- chemid_db$into / (sqrt(sum(chemid_db$into * chemid_db$into)))
    
    dat <- base::merge(pcs_ds,
      chemid_db,
      by = intersect(names(pcs_ds), names(chemid_db)),
      all = TRUE
    )

    ## extract cosine value
    cos_mat <- object@mat
    cos <- cos_mat[match(chemid, rownames(cos_mat)), match(pcs, colnames(cos_mat))]
    cos <- round(cos, digits = 3)

    ## make summary annotation table (rows - unique DB features)
    anno <- object@anno
    anno_mat <- anno[which(anno$ds_peakid %in% pcs_ds$peakid &
      anno$db_peakid %in% chemid_db$peakid), ]

    if (nrow(anno_mat) == 0) {
      warning("chemid and pcs have no matching features.")
      ans <- 0
      while (ans < 1) {
        ans <- readline(
          paste0(
            "Do you wish to proceed with plot generation? Enter Y/N "
          )
        )
        ## catch if input is N/n
        ans <- ifelse((grepl("N", ans) | grepl("n", ans)),
          2, 1
        )
        if (ans == 2) {
          stop("method was stopped.")
        }
      }
    }
    message("pcs intensity values are taken from: ", unique(pcs_ds$filename))
    anno_sum <- base::merge(
      chemid_db[, c("peakid", "mz", "rt", "into_scaled")],
      anno_mat,
      by.x = c("peakid"), by.y = "db_peakid", all = TRUE
    )
    anno_sum <- base::merge(
      setNames(anno_sum, nm = c("DB_peakid", "DB_mz", "DB_rt", "DB_into", "peakid")),
      setNames(pcs_ds[, c("peakid", "mz", "rt", "into_scaled")], nm = c("peakid", "mz", "rt", "into")),
      by.x = c("peakid"), by.y = c("peakid"), all = TRUE
    )
    anno_sum <- anno_sum[, c("DB_mz", "DB_rt", "DB_into", "DB_peakid", "mz", "rt", "into", "peakid")]
    anno_sum <- anno_sum[order(anno_sum$DB_into, decreasing = TRUE), ] ## order by DB features intensity

    ## make spectra mirror plot
    gg_title <- paste0(
      "chemid: ", chemid, " (dbname: ", unique(chemid_db$dbname), ", dbid: ", unique(chemid_db$dbid), ")\n",
      "pcs: ", pcs, " (", unique(pcs_ds$filename), ")\n",
      "cos:", cos
    )
    gg_spectra <- mirrorSPECTRAanno(dat = dat, gg_title = gg_title)

    ## make color_by column
    dat$color_by <- 1
    dat$color_by[is.na(dat$chemid)] <- 2
    gg_labels <- setNames(c(paste("chemid", chemid), paste("PCS", pcs)), nm = c(1:2))
    gg_cols <- c("Black", "#FDE725FF")

    ## make rt deviation plot
    gg_rt <- compareRT(dat = dat, gg_cols = gg_cols, gg_labels = gg_labels, gg_title = gg_title)

    if (!is.null(out_dir)) {
      ## save plot
      ggplot2::ggsave(
        file = paste0("chemid-", chemid, "_pcs-", pcs, "_SPECTRA.png"),
        gg_spectra,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 22, height = 16, units = "cm",
        limitsize = FALSE
      )
      ggplot2::ggsave(
        file = paste0("chemid-", chemid, "_pcs-", pcs, "_RT.png"),
        gg_rt,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 22, height = 16, units = "cm",
        limitsize = FALSE
      )
      write.csv(anno_sum,
        file = paste0(out_dir, "/chemid-", chemid, "_pcs-", pcs, "_summary.csv"),
        row.names = FALSE
      )
    } else {
      gridExtra::grid.arrange(gg_spectra)
      gridExtra::grid.arrange(gg_rt)
    }
    return(anno_sum)
  }
)


# checkADDUCTS ------------------------------------------------------------
## data.frame must contain columns 'mz', 'mzMin', 'mzMax', 'rt', 'rtMin', 'rtMax'
setMethod("checkADDUCTS",
          signature = "massFlowAnno",
          function(object,
                   chemid = NULL,
                   adducts = NULL,
                   pcs = NULL,
                   out_dir = NULL) {
            if (is.null(chemid)) {
              stop("'chemid' is required")
            }
            if (length(chemid) > 1) {
              stop("only one 'chemid' is permitted")
            }
            if (is.null(pcs)) {
              stop("'pcs' is required")
            }
            if (length(pcs) > 1) {
              stop("only one 'pcs' is permitted")
            }
            ## extract features for PCS and CHEMID
            ds <- object@ds
            data_mat <- object@data
            samples_ind <- match(object@samples$filename, colnames(data_mat))
            data_mat <- data_mat[, samples_ind]
            db <- object@db
            
            ## extract features from CHEMID:
            chemid_db <- db[db$chemid == chemid, ]
            chemid_db$into_scaled <- chemid_db$into / (sqrt(sum(chemid_db$into * chemid_db$into)))
            ## find which CHEMID features correspond to validated adducts
            chemid_db$adduct <- findADDUCTS(feats = chemid_db[ , c("mz", "rt")], ## take only numerical columns, otherwise as.matrix() converts to character
                                            adducts = adducts[ , c("mzMin", "mzMax", "rtMin", "rtMax")])
            if (!any(chemid_db$adduct > 0)) {
              return("chemical standards doesn't match any of the provided adducts")
            }
            
            ## extract features for PCS:
            pcs_ds <- ds[ds$pcs == pcs, ]
            pcs_ds$into_scaled <- pcs_ds$into / (sqrt(sum(pcs_ds$into * pcs_ds$into)))
            ## correlate all PCS features with all features
            pcs_cor <- corPCS(pcs_ds = pcs_ds, data_mat = data_mat)
            ## find which PCS features correspond to validated adducts
            pcs_ds$adduct <- findADDUCTS(feats = pcs_ds[ , c("mz", "rt")],
                                         adducts[ , c("mzMin", "mzMax", "rtMin", "rtMax")])
            if (!any(pcs_ds$adduct > 0)) {
              return("pcs doesn't match any of the provided adducts")
            }
            
            ## make summary annotation table (rows - unique DB features)
            anno <- object@anno
            anno_mat <- anno[which(anno$ds_peakid %in% pcs_ds$peakid &
                                     anno$db_peakid %in% chemid_db$peakid), ]
          
            ## extract cosine value
            cos_mat <- object@mat
            cos <- cos_mat[match(chemid, rownames(cos_mat)), match(pcs, colnames(cos_mat))]
            cos <- round(cos, digits = 3)
            
            if (nrow(anno_mat) == 0) {
              return("chemid and pcs have no matching features.")
            }
            message("pcs intensity values are taken from: ", unique(pcs_ds$filename))
            anno_sum <- base::merge(
              chemid_db[, c("peakid", "mz", "rt", "into_scaled", "adduct")],
              anno_mat[ , c("db_peakid", "ds_peakid")],
              by.x = c("peakid"), by.y = "db_peakid", all = TRUE
            )
            anno_sum <- base::merge(
              setNames(anno_sum, nm = c("DB_peakid", "DB_mz", "DB_rt", "DB_into", "adduct", "peakid")),
              setNames(pcs_ds[, c("peakid", "mz", "rt", "into_scaled")], nm = c("peakid", "mz", "rt", "into")),
              by.x = c("peakid"), by.y = c("peakid"), all = TRUE
            )
            anno_sum <- anno_sum[, c("adduct", "DB_mz", "DB_rt", "DB_into", "DB_peakid", "mz", "rt", "into", "peakid")]
            anno_sum <- anno_sum[order(anno_sum$DB_into, decreasing = TRUE), ] ## order by DB features intensity
           
            ## make plots for every adduct
            adduct_inds <- which(pcs_ds$adduct > 0)
            for (ind in adduct_inds) {
              adduct_no <- pcs_ds$adduct[ind]
              if (!adduct_no %in% chemid_db$adduct) {
                next()
              }
              pcs_ds$correlation <- unlist(pcs_cor[[ind]])
              dat <- base::merge(pcs_ds,
                                 chemid_db,
                                 by = intersect(names(pcs_ds), names(chemid_db)),
                                 all = TRUE
              )
              gg_title <- paste0(
                "chemid: ", chemid, " (dbname: ", unique(chemid_db$dbname), ", dbid: ", unique(chemid_db$dbid), ")\n",
                "pcs: ", pcs, " (", unique(pcs_ds$filename), ")\n",
                "cos:", cos, ")\n",
                "adduct: ", adducts$mz[adduct_no]
              )
              gg_spectra <- mirrorSPECTRAtarget(dat = dat, adduct_ind = adduct_no, gg_title = gg_title)

              if (!is.null(out_dir)) {
                ## save plot
                ggplot2::ggsave(
                  file = paste0("chemid-", chemid, "_pcs-", pcs, "_adduct-", adduct_no, "_SPECTRA.png"),
                  gg_spectra,
                  device = "png",
                  path = out_dir,
                  height = 150, width = 200, units = "mm", dpi = "print")
                } else {
                  gridExtra::grid.arrange(gg_spectra)
                }
            }
            write.csv(anno_sum,
                      file = paste0(out_dir, "/chemid-", chemid, "_pcs-", pcs, "_summary.csv"),
                      row.names = FALSE
            )
            return(anno_sum)
          }
)


# comparePCS --------------------------------------------------------------
#' @aliases comparePCS
#'
#' @title Compare selected pseudo chemical spectra.
#'
#' @description Method looks into selected pseudo chemical spectra and returns multiple plots:
#' (1) spectra in the most intense samples; (2) retention time differences; (3) intensity drift during acquisition time.
#' If out_dir is provided, generated plots are written as png files. Otherwise, plots are plotted on current graphical device.
#' \code{massFlowAnno} class object must be annotated using \code{\link{annotateDS}} method first.
#'
#' @param object \code{massFlowAnno} class object.
#' @param pcs \code{numeric} specifying two or more pseudo chemical spectra to look at.
#' @param out_dir \code{character} specifying desired directory for generated figure. Default set to NULL and figure is plotted to graphical device.
#'
#' @return Method returns three plots comparing selected pseudo chemical spectra.
#'
#' @export
#'
setMethod("comparePCS",
  signature = "massFlowAnno",
  function(object,
             pcs = NULL,
             out_dir = NULL) {
    if (is.null(pcs)) {
      stop("'pcs' is required")
    }
    if (length(pcs) == 1) {
      warning("only one pcs was supplied")
    }
    if (length(pcs) > 3) {
      warning("only the first three pcs will be plotted")
      pcs <- pcs[1:3]
    }
    ds <- object@ds
    data_mat <- object@data
    samples_ind <- match(object@samples$filename, colnames(data_mat))
    data_mat <- data_mat[, samples_ind]
    gg_title <- paste0(
      "pcs: ", paste0(pcs, collapse = ", ")
    )
    
    ## for every PCS:
    ## (1) extract intensity values in the most abundant sample
    ## (2) correlate intensities with the main adduct
    pcs_ds <- lapply(pcs, ds = ds, data_mat = data_mat, FUN = prepPCS)
    pcs_ds <- do.call(rbind, pcs_ds)
    
    ## extract intensity matrix by sample run order (samples are already ordered)
    pcs_mat <- lapply(pcs, samples = object@samples, ds = ds, data_mat = data_mat, FUN = extractPCS)
    pcs_mat <- do.call(rbind, pcs_mat)
    pcs_mat$into[pcs_mat$into == 0] <- NA
    
    ## spectra comparison plot
    if (length(pcs) == 2) {
      ## compare one-to-one in a mirror spectra
      gg_spectra <- mirrorSPECTRApcs(dat = pcs_ds, gg_title = gg_title)
    } else {
      ## lsit spectra in separate facets
      gg_spectra <- multipleSPECTRA(dat = pcs_ds, gg_title = gg_title)
    }
    if (length(pcs) == 1) {
      ## color by peakids rather than pcs
      pcs_ds$color_by <- match(pcs_ds$peakid, pcs_ds$peakid)
      gg_labels <- setNames(c(paste("peakid", pcs_ds$peakid)), nm = seq(length(pcs_ds$peakid)))
      gg_cols <- viridis::viridis(length(pcs_ds$peakid))
      pcs_mat$color_by <- match(pcs_mat$peakid, pcs_ds$peakid)
    } else {
      pcs_ds$color_by <- match(pcs_ds$pcs, pcs)
      gg_labels <- setNames(c(paste("PCS", pcs)), nm = seq(length(pcs)))
      gg_cols <- viridis::viridis(length(pcs))
      pcs_mat$color_by <- match(pcs_mat$pcs, pcs)
    }
    gg_rt <- compareRT(dat = pcs_ds, gg_cols = gg_cols, gg_labels = gg_labels, gg_title = gg_title)    
    gg_into <- compareINTENSITY(dat = pcs_mat, gg_cols = gg_cols, gg_labels = gg_labels, gg_title = gg_title)
   
    if (!is.null(out_dir)) {
      ## save plot
      ggplot2::ggsave(
        file = paste0("pcs-", paste0(pcs, collapse = "_"), "_SPECTRA.png"),
        gg_spectra,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 22, height = 16, units = "cm",
        limitsize = FALSE
      )
      ggplot2::ggsave(
        file = paste0("pcs-", paste0(pcs, collapse = "_"), "_RT.png"),
        gg_rt,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 22, height = 16, units = "cm",
        limitsize = FALSE
      )
      ggplot2::ggsave(
        file = paste0("pcs-", paste0(pcs, collapse = "_"), "_INTO.png"),
        gg_into,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 22, height = 16, units = "cm",
        limitsize = FALSE
      )
    } else {
      gridExtra::grid.arrange(gg_spectra)
      gridExtra::grid.arrange(gg_rt)
      gridExtra::grid.arrange(gg_into)
    }
  }
)
