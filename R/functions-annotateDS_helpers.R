getINTENSE <- function(pcs, int_dat, dat) {
  pcs_ind <- which(int_dat$pcs == pcs)
  pcs_sums <- apply(dat[pcs_ind, ], 2, sum)
  sums_max <- base::which.max(pcs_sums)
  data.frame(into = dat[pcs_ind, sums_max],
             filename = names(sums_max),
             stringsAsFactors = FALSE)
  
}


plotSPECTRA <- function(dat, gg_cols, gg_labels, gg_title) {
  gg <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz, y = 0, yend = into, color = as.factor(color_by)),
                          size = 1,  na.rm = TRUE) +
    ggplot2::scale_color_manual(name = "",
                                values = gg_cols, labels = gg_labels) +
    ggplot2::scale_x_continuous(name = "m/z") +
    ggplot2::scale_y_continuous(name = "Intensity") +
    ggplot2::facet_wrap(~color_by, labeller = ggplot2::as_labeller(gg_labels), ncol = 1) + 
    ggplot2::ggtitle(gg_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(size = 0.1)
    ) 
  return(gg)
}

