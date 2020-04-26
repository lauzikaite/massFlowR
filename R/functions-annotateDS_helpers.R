getINTENSE <- function(pcs, int_dat, dat) {
  pcs_ind <- which(int_dat$pcs == pcs)
  pcs_sums <- apply(dat[pcs_ind, ], 2, sum)
  sums_max <- base::which.max(pcs_sums)
  data.frame(
    into = dat[pcs_ind, sums_max],
    filename = names(sums_max),
    stringsAsFactors = FALSE
  )
}

mirrorSPECTRAanno <- function(dat, gg_title) {
  ggplotTHEME(
    gg =
      mirrorSPECTRA(gg = ggplot2::ggplot() +
        ## make top spectra: DB chemid is the same in every facet
        ggplot2::geom_segment(
          data = subset(dat, !is.na(chemid)),
          ggplot2::aes_string(x = "mz", xend = "mz", y = 0, yend = "into_scaled"),
          color = "black",
          size = 1, na.rm = TRUE
        ) +
        ## make the bottom spectra: matching PCS
        ggplot2::geom_segment(
          data = subset(dat, is.na(chemid)),
          ggplot2::aes(x = mz, xend = mz, y = 0, yend = -into_scaled, color = correlation),
          size = 1, na.rm = TRUE
        )),
    gg_title = gg_title
  )
}

mirrorSPECTRAtarget <- function(dat, adduct_ind, gg_title) {
  ggplotTHEME(
    gg =
      mirrorSPECTRA(gg = ggplot2::ggplot() +
                      ## make top spectra: DB chemid
                      ggplot2::geom_segment(
                        data = subset(dat, !is.na(chemid) & adduct != adduct_ind),
                        ggplot2::aes(x = mz, xend = mz, y = 0, yend = into_scaled),
                        color = "black",
                        size = 1, na.rm = TRUE
                      ) +
                      ggplot2::geom_segment(
                        data = subset(dat, !is.na(chemid) & adduct == adduct_ind),
                        ggplot2::aes(x = mz, xend = mz, y = 0, yend = into_scaled),
                        color = "red",
                        size = 1, na.rm = TRUE
                      ) +
                      ## make the bottom spectra: matching PCS
                      ggplot2::geom_segment(
                        data = subset(dat, is.na(chemid)),
                        ggplot2::aes(x = mz, xend = mz, y = 0, yend = -into_scaled, color = correlation),
                        size = 1, na.rm = TRUE
                      )),
    gg_title = gg_title
  )
  
}

mirrorSPECTRApcs <- function(dat, gg_title) {
  ggplotTHEME(
    gg =
      mirrorSPECTRA(gg = ggplot2::ggplot() +
        ## make top spectra: first PCS
        ggplot2::geom_segment(
          data = subset(dat, pcs == unique(dat$pcs)[1]),
          ggplot2::aes(x = mz, xend = mz, y = 0, yend = into_scaled, color = correlation),
          size = 1, na.rm = TRUE
        ) +
        ## make the bottom spectra: second PCS
        ggplot2::geom_segment(
          data = subset(dat, pcs == unique(dat$pcs)[2]),
          ggplot2::aes(x = mz, xend = mz, y = 0, yend = -into_scaled, color = correlation),
          size = 1, na.rm = TRUE
        )),
    gg_title = gg_title
  )
}

mirrorSPECTRA <- function(gg) {
  gg +
    ggplot2::scale_colour_gradientn(
      name = "",
      colours = c("red", "yellow", "green", "lightblue", "darkblue"),
      values = c(1, 0.5, 0, -0.5, -1),
      limits = c(-1, 1)
    ) +
    ggplot2::scale_x_continuous(name = "m/z") +
    ggplot2::scale_y_continuous(
      name = "Intensity, scaled",
      labels = c(1, 0.5, 0, 0.5, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1)
    ) +
    ggplot2::geom_hline(data = data.frame(y = 0), ggplot2::aes(yintercept = y))
}

multipleSPECTRA <- function(dat, gg_title) {
  ggplotTHEME(
    gg =
      mirrorSPECTRA(gg = ggplot2::ggplot(dat) +
        ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz, y = 0, yend = into_scaled, color = correlation),
          size = 1, na.rm = TRUE
        ) +
        facet_wrap(~pcs, ncol = 1)),
    gg_title = gg_title
  )
}

compareRT <- function(dat, gg_cols, gg_labels, gg_title) {
  ggplotTHEME(
    gg = ggplot2::ggplot() +
      ggplot2::geom_point(data = dat,
                          ggplot2::aes(x = rt, y = mz, color = as.factor(color_by))) +
      ggplot2::scale_color_manual(
        name = "",
        values = gg_cols,
        labels = gg_labels
      ) +
      ggplot2::scale_x_continuous(name = "rt") +
      ggplot2::scale_y_continuous(name = "m/z"),
    gg_title = gg_title
  )
}

compareINTENSITY <- function(dat, gg_cols, gg_labels, gg_title) {
  ggplotTHEME(
    gg = ggplot2::ggplot(dat) +
    ggplot2::geom_line(ggplot2::aes(x = run_order, y = log(into), color = as.factor(color_by), group = peakid)) +
    ggplot2::scale_color_manual(
      name = "",
      values = gg_cols,
      labels = gg_labels
    ) +
    ggplot2::scale_x_continuous(name = "Run order") +
    ggplot2::scale_y_continuous(name = "Intensity, log"),
  gg_title = gg_title
  )
}

ggplotTHEME <- function(gg, gg_title) {
  gg +
    ggplot2::ggtitle(gg_title) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(size = 0.1)
    )
}

prepPCS <- function(pc, ds, data_mat) {
  pc_ds <- ds[ds$pcs == pc, ]
  pc_ds$into_scaled <- pc_ds$into / (sqrt(sum(pc_ds$into * pc_ds$into)))

  ## correlate main adduct with all features
  main <- pc_ds$peakid[which.max(pc_ds$into)]
  main_ind <- which(ds$peakid == main)
  other <- pc_ds$peakid[-which.max(pc_ds$into)]
  other_ind <- match(other, ds$peakid)
  other_cor <- lapply(other_ind, function(ind) {
    cor(as.numeric(data_mat[main_ind, ]),
      as.numeric(data_mat[ind, ]),
      use = "complete"
    )
  })
  pc_ds$correlation <- NA
  pc_ds$correlation[match(other, pc_ds$peakid)] <- unlist(other_cor)
  pc_ds$correlation[match(main, pc_ds$peakid)] <- 1
  return(pc_ds)
}

extractPCS <- function(pc, samples, ds, data_mat) {
  ind <- which(ds$pcs == pc)
  peakids <- ds$peakid[ind]
  dat <- lapply(samples$filename, function(sname) {
    data.frame(
      filename = sname,
      run_order = samples$run_order[match(sname, samples$filename)],
      peakid = peakids,
      into = data_mat[ind, match(sname, samples$filename)],
      pcs = pc,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, dat)
}

# findADDUCT <- function(adduct, ds) {
#   inds <- which(ds$mz >= adduct["mzMin"] &
#           ds$mz <= adduct["mzMax"] &
#           ds$rt >= adduct["rtMin"] &
#           ds$rt <= adduct["rtMax"])
#   ds$peakid[inds]
# }

findADDUCTS <- function(feats, adducts) {
  apply(feats, 1, function(feat) {
    matches <- apply(adducts, 1, function(adduct) {
      all(feat["mz"]  >= adduct["mzMin"] &
            feat["mz"]  <= adduct["mzMax"] &
            feat["rt"] >= adduct["rtMin"] &
            feat["rt"] <= adduct["rtMax"])
    })
    if (any(matches)) {
      which(matches)
    } else {
      0
    }
  })
}

corPCS <- function(pcs_ds, data_mat) {
  peakids <- pcs_ds$peakid
  # correlate all features with all features
  lapply(peakids, function(main) {
    main_ind <- which(ds$peakid == main)
    lapply(peakids, function(other) {
      other_ind <- which(ds$peakid == other)
      cor(as.numeric(data_mat[main_ind, ]),
          as.numeric(data_mat[other_ind, ]),
          use = "pairwise.complete.obs"
      )
    })
  })
}
