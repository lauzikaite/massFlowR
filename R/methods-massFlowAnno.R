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
    cos_mat <- do_alignPEAKS(
      ds = ds,
      tmp = db,
      ds_var_name = "pcs",
      tmp_var_name = "chemid",
      mz_err = mz_err,
      rt_err = rt_err,
      bins = 0.1,
      ncores = ncores,
      anno = TRUE
    )

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
        dbname = paste0(dbnames, collapse = ", "),
        similarity = paste0(cos, collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
    ds_to_db <- do.call("rbind", ds_to_db)
    write.csv(ds_to_db, file = file.path(out_dir, "DSannotated.csv"), row.names = FALSE)

    object@db <- db
    object@mat <- cos_mat
    object@anno <- ds_to_db
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

# plotPCS -----------------------------------------------------------------
#' @aliases plotPCS
#' 
#' @title Plot selected pseudo chemical spectra
#' 
#' @description Method returns a plot with the selected single pseudo chemical spectra in the sample in which it was most intense.
#' If anno set to TRUE, spectra of top five annotated database compounds are also plotted.
#' If out_dir is provided, generated plot is written as a png file. Otherwise, spectra is plotted on current graphical device.
#' \code{massFlowAnno} class object must be annotated using \code{\link{annotateDS}} method first.
#' 
#' @param object \code{massFlowAnno} class object.
#' @param pcs \code{numeric} specifying one or more pseudo chemical spectra to look at.
#' @param anno \code{logical} whether spectra of annotated database compounds should be plotted.
#' @param cutoff \code{numeric} specifying spectral similarity score value.
#' @param out_dir \code{character} specifying desired directory for generated figure. Default set to NULL and figure is plotted to graphical device.
#'
#' @return Method returns a plot with the spectra of single pseudo chemical spectra.
#' 
#' @export
#'
setMethod("plotPCS",
  signature = "massFlowAnno",
  function(object,
             pcs = NULL,
             anno = FALSE,
             cutoff = 0,
             out_dir = NULL) {
    cos_mat <- object@mat
    samples <- object@samples
    ds_dat_cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "peakid", "pcs", "into")

    ds <- object@ds
    pcs_dat <- ds[ds$pcs == pcs, ]
    gg_title <- paste0("pseudo chemical spectra: ", pcs)

    ## extract chemids matching to specified pcs
    pcs_mat <- cos_mat[, match(pcs, colnames(cos_mat))]
    pcs_chemid <- names(which(pcs_mat > cutoff))

    ## order matching chemids by their cosines
    db <- object@db
    cos <- cos_mat[pcs_chemid, as.character(pcs)]
    pcs_chemid <- pcs_chemid[order(cos, decreasing = TRUE)]
    cos <- round(cos[order(cos, decreasing = TRUE)], digits = 3)

    ## make summary annotation table
    pcs_db <- db[db$chemid %in% pcs_chemid, ]
    pcs_db_sum <- pcs_db[match(pcs_chemid, pcs_db$chemid), c("chemid", "dbid", "dbname")]
    pcs_db_sum$"spectral similarity" <- cos
    if (nrow(pcs_db_sum) == 0) {
      pcs_db_sum <- data.frame(chemid = NA, dbid = NA, dbname = NA, "spectral similarity" = NA)
    } else {
      if (nrow(pcs_db_sum) > 5) {
        pcs_db_sum <- pcs_db_sum[1:5, ]
      }
    }
    gg_tb <- gridExtra::tableGrob(pcs_db_sum, rows = NULL, theme = gridExtra::ttheme_default(base_size = 8))

    if (anno == TRUE) {
      if (length(pcs_chemid) == 0) {
        stop("PCS have 0 matching compounds with spectral similarity value higher than the selected cutoff")
      } else {
        if (length(pcs_chemid) > 5) {
          stop("PCS compound have > 5 matching compounds with spectral similarity value higher than the selected cutoff.
               To aid visualisation, select higher cutoff value")
        }
      }
      #### ---- Plot spectral comparison
      dat <- merge(pcs_dat, pcs_db, by = intersect(names(pcs_dat), names(pcs_db)), all = TRUE)

      ## make color_by column
      dat$color_by <- match(dat$chemid, pcs_chemid)
      dat$color_by[is.na(dat$color_by)] <- 0

      ## set column for colors: Dataset PCS is 0, Database chemids are numbered from 1 by their cosine
      gg_labels <- setNames(c("PCS", paste0("chemid ", pcs_chemid)), nm = seq(0, length(cos)))
      gg_cols <- viridis::viridis(n = length(pcs_chemid), begin = 0.2, end = 0.95)
      gg_cols <- c("Black", gg_cols)
    } else {
      #### ---- PCS spectra alone
      dat <- pcs_dat
      dat$color_by <- 0
      gg_labels <- setNames("PCS", nm = 0)
      gg_cols <- "Black"
    }

    gg <- plotSPECTRA(dat, gg_cols, gg_labels, gg_title)
    gg <- gridExtra::arrangeGrob(gg, gridExtra::arrangeGrob(gg_tb), ncol = 1, heights = c(4, 1))

    if (!is.null(out_dir)) {
      ## save plot
      ggplot2::ggsave(
        file = paste0("/plotPCS-", pcs, ".png"),
        gg,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 25, height = 25, units = "cm",
        limitsize = FALSE
      )
    } else {
      gridExtra::grid.arrange(gg)
    }
  }
)

# plotCHEMID --------------------------------------------------------------
#' @aliases plotCHEMID
#' 
#' @title Plot selected database chemical compound.
#' 
#' @description Method returns a plot with the selected chemical compound's spectra and top five annotated pseudo chemical spectra in the dataset.
#' If out_dir is provided, generated plot is written as a png file. Otherwise, spectra is plotted on current graphical device.
#' \code{massFlowAnno} class object must be annotated using \code{\link{annotateDS}} method first.
#' 
#' @param object \code{massFlowAnno} class object.
#' @param chemid \code{numeric} specifying chemical compound id to look at.
#' @param cutoff \code{numeric} specifying spectral similarity score value.
#' @param out_dir \code{character} specifying desired directory for output.
#'
#' @return Method returns a plot with the spectra of selected chemical compound and top five annotated pseudo chemical spectra from the dataset.
#' 
#' @export
#'
setMethod("plotCHEMID",
  signature = "massFlowAnno",
  function(object,
             chemid = NULL,
             cutoff = 0.5,
             out_dir = NULL) {
    cos_mat <- object@mat
    int_dat <- object@data
    samples <- object@samples
    ds_dat_cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "peakid", "pcs", "into")

    ## extract PCS matching to specified CHEMID
    chemid_mat <- cos_mat[match(chemid, rownames(cos_mat)), ]
    chemid_pcs <- names(which(chemid_mat > cutoff))

    if (length(chemid_pcs) == 0) {
      stop("chemical compound have 0 matching PCS with spectral similarity value higher than the selected cutoff")
    } else {
      if (length(chemid_pcs) > 5) {
        stop("chemical compound have > 5 matching PCS with spectral similarity value higher than the selected cutoff.
To aid visualisation, select higher cutoff value")
      }
    }

    #### ---- Spectral comparison
    ## plot how the most intense samples match up to this chemical
    ds <- object@ds
    pcs_ds <- ds[ds$pcs %in% chemid_pcs, ]
    db <- object@db
    chemid_db <- db[db$chemid == chemid, ]

    dat <- base::merge(pcs_ds, chemid_db, by = intersect(names(pcs_ds), names(chemid_db)), all = TRUE)

    ## order matching pcs by their cosines
    cos <- cos_mat[as.character(chemid), chemid_pcs]
    chemid_pcs <- chemid_pcs[order(cos, decreasing = TRUE)]
    cos <- round(cos[order(cos, decreasing = TRUE)], digits = 3)

    ## make color_by column
    dat$color_by <- match(dat$pcs, chemid_pcs)
    dat$color_by[is.na(dat$color_by)] <- 0

    ## make summary annotation table
    pcs_sum <- pcs_ds[match(chemid_pcs, pcs_ds$pcs), c("pcs", "filename")]
    pcs_sum$"spectral similarity" <- cos
    gg_tb <- gridExtra::tableGrob(pcs_sum, rows = NULL, theme = gridExtra::ttheme_default(base_size = 8))

    ## set column for colors: Database chemid is 0, Dataset PCS are numbered from 1 by their cosine
    gg_labels <- setNames(c("chemid", paste0("PCS ", chemid_pcs)), nm = seq(0, length(cos)))
    gg_cols <- c("Black", viridis::viridis(n = length(chemid_pcs), begin = 0.2, end = 0.95))
    gg_title <- paste0("chemid: ", chemid, ", dbname: ", unique(chemid_db$dbname))

    gg <- plotSPECTRA(dat, gg_cols, gg_labels, gg_title)
    gg <- gridExtra::arrangeGrob(gg, gridExtra::arrangeGrob(gg_tb), ncol = 1, heights = c(4, 1))

    if (!is.null(out_dir)) {
      ## save plot
      ggplot2::ggsave(
        file = paste0("/plotCHEMID-", chemid, ".png"),
        gg,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 25, height = 25, units = "cm",
        limitsize = FALSE
      )
    } else {
      gridExtra::grid.arrange(gg)
    }
  }
)

# comparePCS --------------------------------------------------------------
#' @aliases comparePCS
#' 
#' @title Plot selected pseudo chemical spectra intensities over acquisition order.
#' 
#' @description Method returns a plot with the selected single pseudo chemical spectra intensities in all samples.
#' If anno set to TRUE, spectra of top five annotated database compounds are also plotted.
#' If out_dir is provided, generated plot is written as a png file. Otherwise, spectra is plotted on current graphical device.
#' \code{massFlowAnno} class object must be annotated using \code{\link{annotateDS}} method first.
#' 
#' @param object \code{massFlowAnno} class object.
#' @param pcs \code{numeric} specifying one or more pseudo chemical spectra to look at.
#' @param out_dir \code{character} specifying desired directory for generated figure. Default set to NULL and figure is plotted to graphical device. 
#'
#' @return Method returns a plot with the intensities of single pseudo chemical spectra.
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
    cos_mat <- object@mat
    int_dat <- object@data
    samples <- object@samples
    ds_dat_cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "peakid", "pcs", "into")

    #### ---- Intensity over time
    pcs_int <- int_dat[int_dat$pcs %in% pcs, ]

    pcs_dat <- lapply(samples$filename, function(sname) {
      data.frame(
        into = pcs_int[, match(sname, colnames(pcs_int))],
        filename = sname,
        run_order = samples$run_order[match(sname, samples$filename)],
        stringsAsFactors = FALSE
      )
    })
    pcs_dat <- do.call("rbind", pcs_dat)
    pcs_dat$peakid <- pcs_int$peakid
    pcs_dat$pcs <- pcs_int$pcs
    ## replace 0s with NAs
    pcs_dat[pcs_dat == 0] <- NA

    gg_cols <- viridis::viridis(n = length(pcs), begin = 0.2, end = 0.95)
    gg_labels <- setNames(paste0("PCS ", pcs), nm = pcs)
    gg_title <- "selected pseudo chemical spectra"

    gg <- ggplot2::ggplot(pcs_dat) +
      ggplot2::geom_line(ggplot2::aes(x = run_order, y = log(into), color = as.factor(pcs), group = peakid)) +
      ggplot2::scale_color_manual(
        name = "",
        values = gg_cols, labels = gg_labels
      ) +
      ggplot2::ggtitle(gg_title) +
      ggplot2::scale_x_continuous(name = "Run order") +
      ggplot2::scale_y_continuous(name = "Intensity, log") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(size = 0.1)
      )

    if (!is.null(out_dir)) {
      ## save plot
      ggplot2::ggsave(
        file = paste0("/comparePCS-", paste0(pcs, collapse = "_"), ".png"),
        gg,
        device = "png",
        path = out_dir,
        dpi = 300,
        width = 25, height = 20, units = "cm",
        limitsize = FALSE
      )
    } else {
      gridExtra::grid.arrange(gg)
    }
  }
)
