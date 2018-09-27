#' Function to find and plot peaks corresponding to spiked reference compounds
#'
#' @param fname A \code{character} specifying filename of the spiked sample.
#' @param pks A \code{DataFrame} class peak table, created for the spiked sample by \emph{pickPEAKS} function.
#' @param raw A \code{OnDiskMSnExp} class object for the spiked sample's spectrum, created by the \emph{MSnbase::readMSData} function.
#' @param spk A \code{DataFrame} containing information about spiked reference compounds, each row is a unique compound.
#' @param out_dir A \code{character} specifying directory where output data will be saved.
#' @param cwt A \code{CentWaveParam} class object with \emph{centWave} parameters that were used to peak-pick the sample of interest. Object can be created by the \emph{xcms::CentWaveParam} function.
#' @param add_xcmsparams A \code{logical} whether \emph{centWave} parameters should be added to the plot. Default set to TRUE.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 3 (sec).
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.005.
#'
#' @return Function returns a \code{DataFrame} with detected peaks for every spiked compound, also plots spectrum for each peak.
#' @export
#'
checkSPIKED <- function(fname, pks, raw, spk, out_dir, cwt, add_xcmsparams = T, rt_err = 3, mz_err = 0.005) {

  ## extract integration values for the datafile of interest
  spk_vals <- spk[,c(match(fname,  colnames(spk)))]

  ## for every compound in the table
  match_df <- data.frame()

  for(i in 1:nrow(spk)) {

    ### extract EIC for the compound based on targetlynx values
    mr <- data.frame(mzmin = as.numeric(spk[i,"mz"]) - mz_err,
                     mzmax = as.numeric(spk[i,"mz"]) + mz_err,
                     rtmin = as.numeric(spk[i,"rt_min"])*60 - rt_err,
                     rtmax = as.numeric(spk[i,"rt_max"])*60 + rt_err,
                     peakid = "spiked", # color by peakid
                     pair = 1, # facet_grid by pair nmrumber
                     stringsAsFactors = F)

    spec <-  extractSPECTRUM(co = mr, raw = raw)

    ## find corresponding peak in the peak table
    pks_match <- pks %>%
      filter(between(rt, mr$rtmin, mr$rtmax)) %>%
      filter(between(mz, mr$mzmin, mr$mzmax)) %>%
      mutate(pair = 1)

    ###--- if single match in the peak table was found
    if (nrow(pks_match) == 1) {

      message("Single match for compound: ", spk[i,"formula"], ", i = ", i)

      ## extract EIC for the corresponding peak using peak's mz and rt range
      spec_pks <- extractSPECTRUM(co = pks_match, raw = raw)
      spec <- rbind(spec, spec_pks)
      gtitle <- paste0(spk[i,"name"], ": peak was picked")

    } else {

      if (nrow(pks_match) == 0) {

        message("No match for compound: ", spk[i,"formula"], ", i = ", i)

        pks_match <- data.frame(peakid = NA, mz = NA, rt = NA, into = NA, stringsAsFactors = F)
        gtitle <- paste0(spk[i,"name"], ": peak was NOT picked")

      } else {

        message("More than one match for compound: ", spk[i,"formula"], ", i = ", i)

        ## if more than one match, take the peak closest in rt
        rt_close <- which(abs(pks_match$rt - spk[i,"rt"])  == min(abs(pks_match$rt - spk[i,"rt"])))
        pks_match <- pks_match[rt_close,]
        spec_pks <- extractSPECTRUM(co = pks_match, raw = raw)
        spec <- rbind(spec, spec_pks)
        gtitle <- paste0(spk[i,"name"], ": more than one peak picked, closest match")

      }
    }

    ###---- universal part for plotting and output
    ## named vector for set color generation, so that "spiked" is always purple, despite how many peaks
    cols <- setNames(viridis::viridis(end = 0.8, alpha = 0.8,
                                      n = length(unique(spec$peakid))),
                     nm = c("spiked", spec %>% filter(peakid != "spiked") %>% distinct(peakid) %>% pull()))

    g <- ggplot2::ggplot(data = spec) +
      ggplot2::geom_segment(aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = peakid), na.rm = TRUE) +
      ggplot2::geom_point(aes(x = rtime, y = intensity, colour = peakid), na.rm = TRUE) +
      ggplot2:: scale_color_manual(values = cols, name = "") +
      ggplot2::ylab("Intensity") +
      ggplot2::xlab("Retention time") +
      ggplot2::ggtitle(gtitle) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
      ggplot2::theme(legend.position = "bottom")

    if (add_xcmsparams == T & !missing(cwt)) {
        g <- g +
          ggplot2::geom_hline(aes(yintercept = cwt@prefilter[2],
                         linetype = paste0("Prefilter: c(", cwt@prefilter[1], "," , cwt@prefilter[2], ")")),
                     colour= "#1b7837") +
          ggplot2::geom_hline(aes(yintercept = cwt@noise,
                         linetype = paste("Noise:", cwt@noise)),
                     colour= '#762a83') +
          ggplot2::scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                                guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))
        }

    pdf(width = 8, height = 8,
        file = paste0(out_dir, "/", fname, "_spikedMR-", spk[i,"formula"], ".pdf"))
    grid::grid.draw(rbind(ggplotGrob(g), size = "last"))
    dev.off()

    ## save output
    match <- data.frame(sample_id = fname,
                        spiked_id = spk[i,"id"],
                        # spiked_molw = spk[i,"mol_weight"],
                        spiked_mz = spk[i,"mz"],
                        spiked_rt = spk[i,"rt"],
                        pks_id = pks_match$peakid,
                        pks_mz = pks_match$mz,
                        pks_rt = pks_match$rt,
                        pks_into = pks_match$into,
                        pks_into_log = log(pks_match$into),
                        spiked_val = spk_vals[i],
                        spiked_val_log = log(spk_vals[i]))


    match_df <- rbind(match_df, match)
  }
  write.csv(match_df, paste0(out_dir, "/", fname, "_checkSPIKED_output.csv"), quote = F, row.names = F)
  return(match_df)
}


## check list pf pre-selected peaks (can be all from one COMP, can be simply co-eluting, can be from same CLS)
## build correlation  network
## for selected pairs of the network, plot spectra and output cor values and MZ differences
checkPEAKS <- function(pks_mat, pks, spiked, spiked_p, eic, raw, out_dir, cwt, add_xcmsparams = F, thr = 0.95, use_thr = T) {

  ## get cor values between all co-eluting peaks, pair-wise
  peakid_all <- pks_mat  %>%
    pull(peakid)

  pks_cor <- buildCOR(co_ind = peakid_all, eic = eic, pearson = T)

  ## for every peak-pair, extract cor and mz_diff
  pks_cor_mz <- pks_cor %>%
    group_by(x, y) %>%
    ## remove pair duplication
    mutate(x_y = paste(min(c(x, y)), max(c(x, y)), sep = "_")) %>%
    group_by(x_y) %>%
    slice(1) %>%
    mutate(mz_x = pks_mat %>% filter(peakid == x) %>% pull(mz),
           mz_y = pks_mat %>% filter(peakid == y) %>% pull(mz)) %>%
    mutate(mz_dif = round(abs(mz_x - mz_y), digits = 2)) %>%
    ungroup()

  ## extract EIC values for every peakid
  spec <- pks_mat %>%
    group_by(peakid) %>%
    do(extractSPECTRUM(co = ., raw = raw)) %>%
    ungroup()

  ## plot all peaks in one
  ## labels will be mz, color order also mz
  spec <- spec %>%
    rename(spec_peakid = peakid) %>%
    group_by(spec_peakid) %>%
    mutate(mz = pks_mat %>% filter(peakid %in% spec_peakid) %>% distinct(mz) %>% pull() %>% round(., digits = 3)) %>%
    arrange(mz)

  colors <- setNames(viridis::viridis(end = 0.9, ## avoid bright yellow
                                      begin = 0.2, ## avoid dark purple
                                      n = length(unique(spec$mz))),
                              nm = spec %>% distinct(mz) %>% pull())

  g <- ggplot2::ggplot(data = spec) +
    ggplot2::geom_segment(aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(mz)), na.rm = TRUE) +
    ggplot2::geom_point(aes(x = rtime, y = intensity, colour = as.factor(mz)), na.rm = TRUE) +
    ggplot2::scale_color_manual(values = colors,
                       name = "Peak mz") +
    ggplot2::ylab("Intensity") +
    ggplot2::xlab("Retention time") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))

  grid::grid.newpage()
  pdf(width = 8, height = 8,
      file = paste0(out_dir, "/", fname, "_peaks-",paste0(pks_mat$peakid, collapse = "-"),  ".pdf"))
  grid::grid.draw(rbind(ggplotGrob(g), size = "last"))
  dev.off()

  buildNETWORK(poi_co_cor = pks_cor,
               peakgroups = pks_mat,
               co_ind = peakid_all,
               thr = thr,
               pks = pks,
               p = spiked_p,
               plot = T,
               out_dir = out_dir,
               fname = fname,
               return = F,
               colors = colors)

 return(pks_cor_mz)
}

plotPEAKpair <- function(pair, spiked, spec, add_xcmsparams, cwt, out_dir, fname){

  ## extract pair of interest
  spec <- spec %>%
    filter(peakid == pair$x | peakid == pair$y)

  if (spiked == T) {
    spec <- spec %>%
      mutate(peakid = ifelse(peakid == spiked_p, "spiked", peakid))
  }

  cols <- setNames(viridis::viridis(end = 0.8, alpha = 0.8,
                                    n = length(unique(spec$peakid))),
                   nm = if (spiked == T & "spiked" %in% spec$peakid) {
                     c("spiked", spec %>% filter(peakid != "spiked") %>% distinct(peakid) %>% pull()) } else {
                       c(spec %>% filter(peakid != "spiked") %>% distinct(peakid) %>% pull()) }
  )

  g <- ggplot2::ggplot(data = spec) +
    ggplot2::geom_segment(aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(peakid)), na.rm = TRUE) +
    ggplot2::geom_point(aes(x = rtime, y = intensity, colour = as.factor(peakid)), na.rm = TRUE) +
    ggplot2::scale_color_manual(values = cols, name = "") +
    ggplot2::ylab("Intensity") +
    ggplot2::xlab("Retention time") +
    ggplot2::ggtitle(pair$gtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = "bottom")

  if (add_xcmsparams == T & !missing(cwt)) {
    g <- g +
      ggplot2::geom_hline(aes(yintercept = cwt@prefilter[2],
                     linetype = paste0("Prefilter: c(", cwt@prefilter[1], "," , cwt@prefilter[2], ")")),
                 colour= "#1b7837") +
      ggplot2::geom_hline(aes(yintercept = cwt@noise,
                     linetype = paste("Noise:", cwt@noise)),
                 colour= '#762a83') +
      ggplot2::scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                            guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))
  }

  grid::grid.newpage()
  pdf(width = 10, height = 10,
      file = paste0(out_dir, "/", fname, "_peakPair-",pair$x, "-", pair$y, ".pdf"))
  grid::grid.draw(rbind(ggplotGrob(g), size = "last"))
  dev.off()

  return(spec)
}

#' Extract EIC for a chromatographic peak-of-interest using its mz and rt values
#'
#' @param co A \code{DataFrame} containing rtmin, rtmax, mzmin and mzmax values for peak EIC extraction. Must also include column 'peakid'.
#' @param raw A \code{OnDiskMSnExp} class object for the spectrum-of-interest.
#'
#' @return Function returns \code{DataFrame} containing intensity value for each scan in the desired mz/rt range.
#' @export
#'
#' @examples
extractSPECTRUM <- function(co, raw) {

  eic_co <- xcms::chromatogram(raw, aggregationFun = "sum",
                               rt = matrix(c(co$rtmin, co$rtmax), ncol = 2),
                               mz = matrix(c(co$mzmin, co$mzmax), ncol = 2))

  spec <- data.frame(rtime = as.numeric(xcms::rtime(eic_co[1])),
             intensity = as.numeric(xcms::intensity(eic_co[1]))) %>%
    mutate(peakid = co$peakid)

  return(spec)

}

