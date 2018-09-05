#' Function to find and plot peaks corresponding to spiked reference compounds
#'
#' @param fname A \code{character} specifying filename of the spiked sample.
#' @param pks A \code{DataFrame} class peak table, created for the spiked sample by \emph{pickPEAKS} function.
#' @param raw A \code{OnDiskMSnExp} class object for the spiked sample's spectrum, created by the \emph{MSnbase::readMSData} function.
#' @param spk A \code{DataFrame} containing information about spiked reference compounds, each row is a unique compound.
#' @param out_dir A \code{character} specifying directory where output data will be saved.
#' @param paramCWT A \code{CentWaveParam} class object with \emph{centWave} parameters that were used to peak-pick the sample of interest. Object can be created by the \emph{xcms::CentWaveParam} function.
#' @param add_xcmsparams A \code{logical} whether \emph{centWave} parameters should be added to the plot.
#' @param rt_err A \code{numeric} specifying the window for peak matching in the RT dimension. Default set to 3 (sec).
#' @param mz_err A \code{numeric} specifying the window for peak matching in the MZ dimension. Default set to 0.005.
#'
#' @return Function returns a \code{DataFrame} with detected peaks for every spiked compound, also plots spectrum for each peak.
#' @export
#'
#' @examples
checkSPIKED <- function(fname, pks, raw, spk, out_dir, paramCWT, add_xcmsparams = T, rt_err = 3, mz_err = 0.005) {

  ## extract integration values for the datafile of interest
  spk_vals <- spk[,c(match(fname,  colnames(spk)))]

  ## for every compound in the table
  match_df <- data.frame()

  for(i in 1:nrow(spk)) {

    ### extract EIC for the MR compound based on targetlynx values
    mr <- data.frame(mzmin = as.numeric(spk[i,"mz"]) - mz_err,
                     mzmax = as.numeric(spk[i,"mz"]) + mz_err,
                     rtmin = as.numeric(spk[i,"rt_min"])*60 - rt_err,
                     rtmax = as.numeric(spk[i,"rt_max"])*60 + rt_err,
                     pno = "spiked", # color by pno
                     pair = 1, # facet_grid by pair nmrumber
                     stringsAsFactors = F)

    spec <-  massFlowR::extractSPECTRUM(co = mr, raw = raw)

    ## find corresponding peak in the peak table
    pks_match <- pks %>%
      filter(between(rt, mr$rtmin, mr$rtmax)) %>%
      filter(between(mz, mr$mzmin, mr$mzmax)) %>%
      mutate(pair = 1) %>%
      rename(pno = dplyr::matches("order_ID|pid|pno")) ## in case old version peak tables are to be used

    ###--- if single match in the peak table was found
    if (nrow(pks_match) == 1) {

      message("Single match for compound: ", spk[i,"formula"], ", i = ", i)

      spec_pks <- massFlowR::extractSPECTRUM(co = pks_match, raw = raw)
      spec <- rbind(spec, spec_pks)
      gtitle <- paste0(spk[i,"name"], ": peak was picked")

    } else {

      if (nrow(pks_match) == 0) {

        message("No match for compound: ", spk[i,"formula"], ", i = ", i)

        pks_match <- data.frame(pno = NA, mz = NA, rt = NA, into = NA, stringsAsFactors = F)
        gtitle <- paste0(spk[i,"name"], ": peak was NOT picked")

      } else {

        message("More than one match for compound: ", spk[i,"formula"], ", i = ", i)

        ## if more than one match, take the peak closest in rt
        rt_close <- which(abs(pks_match$rt - spk[i,"rt"])  == min(abs(pks_match$rt - spk[i,"rt"])))
        pks_match <- pks_match[rt_close,]
        spec_pks <- massFlowR::extractSPECTRUM(co = pks_match, raw = raw)
        spec <- rbind(spec, spec_pks)
        gtitle <- paste0(spk[i,"name"], ": more than one peak picked, closest match")

      }
    }

    ###---- universal part for plotting and output
    ## named vector for set color generation, so that "spiked" is always purple, despite how many peaks
    cols <- setNames(viridis::viridis(end = 0.8, alpha = 0.8,
                                      n = length(unique(spec$pno))),
                     nm = c("spiked", spec %>% filter(pno != "spiked") %>% distinct(pno) %>% pull()))

    g <- ggplot(data = spec) +
      geom_segment(aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = pno), na.rm = TRUE) +
      geom_point(aes(x = rtime, y = intensity, colour = pno), na.rm = TRUE) +
      scale_color_manual(values = cols, name = "") +
      ylab("Intensity") +
      xlab("Retention time") +
      ggtitle(gtitle) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = "bottom")

    if (add_xcmsparams == T & !missing(paramCWT)) {
        g <- g +
          geom_hline(aes(yintercept = paramCWT@prefilter[2],
                         linetype = paste0("Prefilter: c(", paramCWT@prefilter[1], "," , paramCWT@prefilter[2], ")")),
                     colour= "#1b7837") +
          geom_hline(aes(yintercept = paramCWT@noise,
                         linetype = paste("Noise:", paramCWT@noise)),
                     colour= '#762a83') +
          scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                                guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))
        }

    pdf(width = 6, height = 6,
        file = paste0(out_dir, "/", fname, "_spikedMR-", spk[i,"formula"], ".pdf"))
    grid::grid.draw(rbind(ggplotGrob(g), size = "last"))
    dev.off()

    ## save output
    match <- data.frame(sample_id = fname,
                        spiked_id = spk[i,"id"],
                        spiked_molw = spk[i,"mol_weight"],
                        spiked_mz = spk[i,"mz"],
                        spiked_rt = spk[i,"rt"],
                        pks_id = pks_match$pno,
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
## build correlatio  network
## for selected pairs of the network, plot spectra and output cor values and MZ differences
checkPEAKS <- function(pks_mat, pno_main, eic, raw, out_dir, paramCWT, add_xcmsparams = F, thr = 0.95, use_thr = T) {

  ## get cor values between all co-eluting peaks, pair-wise
  pks_mat <- pks_mat %>%
    rename(pno = dplyr::matches("order_ID|pid|pno"))

  pno_all <- pks_mat  %>%
    pull(pno)

  pks_cor <- buildCOR(co_ind = pno_all, eic = eic, pearson = T)

  ## if only peaks above cor threshold is desired
  if (use_thr == T) {
    pks_cor  <- pks_cor %>%
      filter(cor > thr)
    ## update peak numbers
    pno_all <- pks_cor %>%
      select(x,y) %>%
      pull() %>%
      unique()

    ## update peak table
    pks_mat <- pks_mat %>%
      filter(pno %in% pno_all)
  }

  ## for every peak-pair, extract cor and mz_diff
  pks_cor_mz <- pks_cor %>%
    group_by(x, y) %>%
    ## remove pair duplication
    mutate(x_y = paste(min(c(x, y)), max(c(x, y)), sep = "_")) %>%
    group_by(x_y) %>%
    slice(1) %>%
    mutate(mz_x = pks_mat %>% filter(pno == x) %>% pull(mz),
           mz_y = pks_mat %>% filter(pno == y) %>% pull(mz)) %>%
    mutate(mz_dif = round(abs(mz_x - mz_y), digits = 2)) %>%
    ungroup()

  ## extract EIC values for every pno
  spec <- pks_mat %>%
    group_by(pno) %>%
    do(extractSPECTRUM(co = ., raw = raw))

  g <- pks_cor_mz %>%
    group_by(x,y) %>%
    filter(x == 4, y ==1) %>%
    mutate(gtitle = paste0("Cor: ", round(cor, digits = 2), ". MZ diff: ", mz_dif)) %>%
    do(plotPEAKpair(pair = ., spec = spec, add_xcmsparams = T, paramCWT = paramCWT))






}

plotPEAKpair <- function(pair, spiked, spec, add_xcmsparams, paramCWT){

  ## extract pair of interest
  spec <- spec %>%
    filter(pno == pair$x | pno == pair$y)

  cols <- setNames(viridis::viridis(end = 0.8, alpha = 0.8,
                                    n = length(unique(spec$pno))),
                   nm = if (spiked == T) {
                     c("spiked", spec %>% filter(pno != "spiked") %>% distinct(pno) %>% pull()) } else {
                       c(spec %>% filter(pno != "spiked") %>% distinct(pno) %>% pull()) }
  )

  g <- ggplot(data = spec) +
    geom_segment(aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(pno)), na.rm = TRUE) +
    geom_point(aes(x = rtime, y = intensity, colour = as.factor(pno)), na.rm = TRUE) +
    scale_color_manual(values = cols, name = "") +
    ylab("Intensity") +
    xlab("Retention time") +
    ggtitle(pair$gtitle) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom")

  if (add_xcmsparams == T & !missing(paramCWT)) {
    g <- g +
      geom_hline(aes(yintercept = paramCWT@prefilter[2],
                     linetype = paste0("Prefilter: c(", paramCWT@prefilter[1], "," , paramCWT@prefilter[2], ")")),
                 colour= "#1b7837") +
      geom_hline(aes(yintercept = paramCWT@noise,
                     linetype = paste("Noise:", paramCWT@noise)),
                 colour= '#762a83') +
      scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                            guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))
  }

  return(g)

}




checkSPIKEDcomps <- function(i, fname, pks, raw, eic, spk_pks, cor_mat, out_dir, paramCWT, rt_err = 3, mz_err = 0.005, thr = 0.95, pearson = TRUE) {

  if (missing(i)) {
    ## for all spiked MR compounds
    i_list <- 1:nrow(spk_pks)

  } else {
    ## for selected spiked MR compound only
    i_list <- i
  }

  pks <- pks %>%
    ## rename pid to pno, if old-version table is being used
    rename(pno = matches("pid|pno"))

  cor_data <- data.frame()
  for(i in i_list) {
    message("Checking spiked compound: ", spk_pks[i,"spiked_id"], ", i = ", i)

    ## get the peak ID which represents this compound
    poi <- spk_pks[i,"pks_id"]
    cor_data_i <- checkANYpeak(p_list = poi, fname = fname, pks = pks, raw = raw, eic = eic, out_dir = out_dir, paramCWT = paramCWT, use_spk = T, use_thr = T, add_xcmsparams = F)

    ## get mz difference to the absolute spiked MZ
    cor_data_i  <- data.frame(
      mz_molw_dif = spk_pks[i,"spiked_molw"] - (cor_data_i %>% distinct(mz_x, mz_y) %>% pull()),
      mz_m_dif = spk_pks[i,"spiked_mz"] - (cor_data_i %>% distinct(mz_x, mz_y) %>% pull()),
      stringsAsFactors = F) %>%
      mutate(spiked_i = i)

    cor_data <- rbind(cor_data, cor_data_i)
  }

  return(cor_data)
}


checkANYpeak <- function(p_list, fname, pks, raw, eic, out_dir, paramCWT, use_spk, use_thr = T, thr = 0.95, add_xcmsparams = F) {

  cor_data <- data.frame()

  for(p in p_list) {

    message("Checking peak with pno: ", p )

    ## get component peaks for the poi
    p_comp <- pks[p, "comp"]
    p_co <- pks %>%
      filter(comp == p_comp) %>%
      mutate(poi = ifelse(pno == p, T, F)) %>%
      ## main peak - the spiked one, or the most intense
      mutate(main =  if (use_spk == T) { ifelse(pno == p, T, F) } else { ifelse(into == max(into), T, F) })

    if (nrow(p_co) == 1) { next("No co-eluting peaks for PID: ", p) }

    ## get cor values between all co-eluting peaks, pair-wise
    p_co_pno <- pull(p_co, pno)
    p_co_cor <- buildCOR(co_ind = p_co_pno, eic = eic, pearson = T)

    ## if only peaks above cor threshold is desired
    if (use_thr == T) {
      p_co_cor <- p_co_cor %>%
        filter(cor > thr)
      p_co <- p_co %>%
        filter(pno %in%
                 (p_co_cor %>%
                 select(x,y) %>%
                 pull()))
      }

    ## extract EIC values for every corelated peak
    spec <- p_co %>%
      group_by(pno) %>%
      do(extractSPECTRUM(co = ., raw = raw)) %>%
      # mutate(co_pno = pno) %>%
      # mutate(poi = poi) %>%
      # mutate(main = main) %>%
      ungroup()

    ##  save MZ differences within the component
    p_co_cor <- p_co_cor %>%
      group_by(x, y) %>%
      mutate(x_y = paste(min(c(x, y)), max(c(x, y)), sep = "_")) %>%
      group_by(x_y) %>%
      slice(1) %>%
      mutate(mz_x = pks %>% filter(pno == x) %>% pull(mz),
             mz_y = pks %>% filter(pno == y) %>% pull(mz)) %>%
      mutate(mz_dif = round(abs(mz_x - mz_y), digits = 2)) %>%
      ungroup()

    ## save POI-co peak pairs
    cor_data <- rbind(cor_data,
                      p_co_cor %>%
                        filter(x == p | y == p))

    ####---- plot EICs of main adduct and all co-eluting peaks
    ## only works well if peaks are only correlated with the POI
    ## need to re-write for more generic cases
    # g_labels <- p_co_cor %>%
    #   filter(x == p | y == p) %>%
    #   # mutate(not_p_pno = ifelse(x == p, y, x)) %>%
    #   # group_by(not_p_pno) %>%
    #   # slice(1) %>%
    #   mutate(g_text = paste0("Cor: ", round(cor, digits = 2), ". MZ diff: ", mz_dif)) %>%
    #   ungroup() %>%
    #   select(x, y, not_p_pno, g_text) %>%
    #   bind_rows(data.frame(x = p, y = p, not_p_pno = p, g_text = paste("Peak-of-Interest"), stringsAsFactors = F)) %>%
    #   arrange(not_p_pno)
    #
    # g <- ggplot(spec, aes(x = rtime, xend = rtime, y = 0, yend = intensity)) +
    #   ## exclusion of column 'co_pid' allows to retain data points in each facet
    #   geom_segment(data = filter(spec, main == T) %>% select(., -co_pno),
    #                alpha = 0.6, colour = "black", na.rm = TRUE) +
    #   geom_point(data = filter(spec, main == T)  %>% select(., -co_pno), aes(y = intensity),
    #              alpha = 0.6, colour = "black", na.rm = TRUE) +
    #   geom_segment(data = filter(spec, main == F),
    #                aes(colour = as.factor(co_pno)),
    #                alpha = 0.8, na.rm = TRUE) +
    #   geom_point(data = filter(spec, main == F),
    #              aes(y = intensity,colour = as.factor(co_pno)),
    #              alpha = 0.8, na.rm = TRUE) +
    #   scale_color_viridis_d(end = 0.8,
    #                         alpha = 0.8,
    #                         name = "Peak ID") +
    #   facet_wrap(~co_pno,
    #              labeller =  as_labeller(setNames(g_labels$g_text, nm = g_labels$not_p_pno)),
    #              scales = "free_y") +
    #   xlab("Retention time") +
    #   ylab("Intensity") +
    #   theme_bw() +
    #   theme(legend.position = "bottom")
    #
    # if (add_xcmsparams == T) {
    #   g <- g +
    #     geom_hline(aes(yintercept = paramCWT@prefilter[2],
    #                    linetype = paste0("Prefilter: c(", paramCWT@prefilter[1], "," , paramCWT@prefilter[2], ")")),
    #                colour= "#1b7837") +
    #     geom_hline(aes(yintercept = paramCWT@noise,
    #                    linetype = paste("Noise:", paramCWT@noise)),
    #                colour= '#762a83') +
    #     scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
    #                           guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))
    #   }
    #
    #   grid::grid.newpage()
    #   pdf(width = 10, height = 10,  paper = "a4",
    #       file = paste0(out_dir, "/", fname, "_poi-", p, "_allCO_spectra.pdf"))
    #   grid::grid.draw(ggplotGrob(g))
    #   dev.off()


  }

  return(cor_data)
}










#' Extract EIC for a peak-of-interest using its range of mz and rt values
#'
#' @param co A \code{DataFrame} containing rtmin, rtmax, mzmin and mzmax values for peak EIC extraction. Must also include columns 'pno' and 'pair'.
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
    mutate(pno = co$pno)
           # pair = co$pair)

  return(spec)

}

