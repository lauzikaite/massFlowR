#' Title
#'
#' @param fname
#' @param pks
#' @param raw
#' @param eic
#' @param spk
#' @param out_dir
#' @param prefilter
#' @param noise
#' @param rt_err
#' @param mz_err
#'
#' @return
#' @export
#'
#' @examples
checkSPIKED <- function(fname, pks, raw, eic, spk, out_dir, prefilter, noise, rt_err = 3, mz_err = 0.005) {

  ####---- for all spiked MR compounds ----
  match_df <- data.frame()

  for(i in 1:nrow(spk)) {

    ## match spiked dataset and perak table for the single datafile (i.e. fname)
    spk_match <- spk[,c(match(fname,  colnames(spk)))]
    spiked_val <- as.numeric(spk_match[i])
    spiked_val_log <- log(as.numeric(spk_match[i]))

    ## get mz/rt windows for matching
    mz <- as.numeric(spk[i,"mz"])
    mzr <- matrix(c(mz - mz_err, mz + mz_err), ncol = 2, byrow = TRUE)

    rt <- as.numeric(spk[i,"rt"])*60
    rt_min <- as.numeric(spk[i,"rt_min"])*60
    rt_max <- as.numeric(spk[i,"rt_max"])*60
    rtr <- matrix(c(rt_min - rt_err, rt_max + rt_err), ncol = 2, byrow = TRUE)

    ### extract EIC for full matching spectrum based on targetlynx values
    eic_i <-  xcms::chromatogram(raw,
                                 aggregationFun = "sum",
                                 rt = rtr,
                                 mz = mzr)
    eic_int <- xcms::intensity(eic_i[1,1])
    eic_rtime <- xcms::rtime(eic_i[1,1])
    gdf <- data.frame(rtime = as.numeric(eic_rtime),
                      intensity = as.numeric(eic_int),
                      peak = rep(1, length(eic_rtime)))

    ## find corresponding peak in the peak table
    pks_match <- pks %>%
      filter(between(rt, rtr[,1], rtr[,2])) %>%
      filter(between(mz, mzr[,1], mzr[,2]))

    pks_match_id <- pks_match %>%
      select(dplyr::matches("order_ID|pid")) %>%
      pull()

    ## extract picked peak intensities
    pks_val <- as.numeric(pks_match$into)
    pks_val_log <- log(as.numeric(pks_match$into))

    ###--- if single match in the peak table ----
    if (nrow(pks_match) == 1) {

      message("Single match: id = ", spk[i,"id"], ", i = ", i)

      eic_pks <- xcms::chromatogram(raw, aggregationFun = "sum",
                                    rt = c(pks_match$rtmin, pks_match$rtmax),
                                    mz = c(pks_match$mzmin, pks_match$mzmax))
      eic_pks_int <- xcms::intensity(eic_pks[1,1])
      eic_pks_rtime <- xcms::rtime(eic_pks[1,1])
      gdf <- rbind(gdf, data.frame(rtime = as.numeric(eic_pks_rtime),
                                  intensity = as.numeric(eic_pks_int),
                                  peak = rep(2, length(eic_pks_rtime))))
      cmp <- paste0(spk[i,"name"], ": peak was picked")

      g <- ggplot(gdf) +
        geom_segment(data = gdf,  aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(peak)), na.rm = TRUE) +
        geom_point(data = gdf, aes(x = rtime, y = intensity, colour = as.factor(peak)), na.rm = TRUE) +
        scale_colour_brewer(palette="Paired",
                          name = "", labels = c("Full spectrum", "Picked peak")) +
        geom_hline(aes(yintercept= prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
        geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
        scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                              guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))

      ## save output
      match <- data.frame(sample_id = fname,
                          spiked_id = spk[i,"id"],
                          spiked_molw = spk[i,"mol_weight"],
                          spiked_mz = mz,
                          spiked_rt = rt,
                          pks_id = pks_match_id,
                          pks_mz = pks_match$mz,
                          pks_rt = pks_match$rt,
                          pks_into = pks_val,
                          pks_into_log = pks_val_log,
                          spiked_val = spiked_val,
                          spiked_val_log = spiked_val_log)



    } else {

      if (nrow(pks_match) == 0) {

        message("No match: id = ", spk[i,"id"], ", i =", i)

        cmp <- paste0(spk[i,"name"], ": peak was NOT picked")
        g <- ggplot(gdf) +
          geom_segment(data = gdf,  aes(x = rtime, xend = rtime, y = 0, yend = intensity), colour = "#a6cee3", na.rm = TRUE) +
          geom_point(data = gdf, aes(x = rtime, y = intensity), colour = "#a6cee3", na.rm = TRUE) +
          geom_hline(aes(yintercept= prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
          geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
          scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                                guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))

        match <- data.frame(sample_id = fname,
                            spiked_id = spk[i,"id"],
                            spiked_molw = spk[i,"mol_weight"],
                            spiked_mz = mz,
                            spiked_rt = rt,
                            pks_id = NA,
                            pks_mz = NA,
                            pks_rt = NA,
                            pks_into = NA,
                            pks_into_log = NA,
                            spiked_val = spiked_val,
                            spiked_val_log = spiked_val_log)

      } else {

        message("More than one match: id = ", spk[i,"id"], ", i = ",i) ## (!) duplicated matches removal to be updated
        rt_close <- which(abs(pks_match$rt - rt)  == min(abs(pks_match$rt - rt)))
        if (fname == "AIRWAVE_LNEG_ToF06_P73W84_SR" & i == 7 ) { rt_close <- 1 } # cheap hack, need to use sample specific rt
        eic_pks <- xcms::chromatogram(raw, aggregationFun = "sum",
                                      rt = c(pks_match$rtmin[rt_close], pks_match$rtmax[rt_close]),
                                      mz = c(pks_match$mzmin[rt_close], pks_match$mzmax[rt_close]))
        eic_pks_int <- xcms::intensity(eic_pks[1,1])
        eic_pks_rtime <- xcms::rtime(eic_pks[1,1])
        gdf <- rbind(gdf, data.frame(rtime = as.numeric(eic_pks_rtime),
                                     intensity = as.numeric(eic_pks_int),
                                     peak = rep(2, length(eic_pks_rtime))))
        cmp <- paste0(spk[i,"name"], ": more than one peak picked, closest match")

        g <- ggplot(gdf) +
          geom_segment(data = gdf,  aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(peak)), na.rm = TRUE) +
          geom_point(data = gdf, aes(x = rtime, y = intensity, colour = as.factor(peak)), na.rm = TRUE) +
          scale_colour_brewer(palette="Paired",
                              name = "", labels = c("Full spectrum", "Picked peak")) +
          geom_hline(aes(yintercept= prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
          geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
          scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                                guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837"))))


        match <- data.frame(sample_id = fname,
                            spiked_id = spk[i,"id"],
                            spiked_molw = spk[i,"mol_weight"],
                            spiked_mz = mz,
                            spiked_rt = rt,
                            pks_id = pks_match_id[rt_close],
                            pks_mz = pks_match$mz[rt_close],
                            pks_rt = pks_match$rt[rt_close],
                            pks_into = pks_val[rt_close],
                            pks_into_log = pks_val_log[rt_close],
                            spiked_val = spiked_val,
                            spiked_val_log = spiked_val_log)
      }
    }

    ####---- finish with plotting ----
    g <- g +
      ylab("Intensity") +
      xlab("Retention time") +
      ggtitle(cmp) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(filename = paste0(out_dir, "/", fname, "_spikedMR-", spk[i,"formula"], ".png"), plot = g,width = 10, height = 8, units = "in")

    match_df <- rbind(match_df, match)
  }
  write.csv(match_df, paste0(out_dir, "/", fname, "_checkSPIKED_output.csv"), quote = F, row.names = F)
  return(match_df)
}


checkSPIKEDcomps <- function(i, fname, pks, raw, eic, spk_pks, cor_mat, out_dir, paramCWT, rt_err = 3, mz_err = 0.005, thr = 0.95) {

  if (missing(i)) {

    ## for all spiked MR compounds
    i_list <- 1:nrow(spk_pks)

  } else {

    ## for selected spiked MR compound only
    i_list <- i
  }

  for(i in i_list) {

    message("Checking spiked compound: ", spk[i,"id"], ", i = ", i)

    ## get all peaks in the same component as the matched peak
    spk_pks_comp <- pks %>%
      filter(comp == (pks %>% filter(pid == spk_pks[i,"pks_id"]) %>% select(comp) %>% pull) ) %>%
      mutate(spiked = ifelse(pid == spk_pks[i,"pks_id"], T, F)) %>%
      ## which peak is the most intense - i.e. the main
      mutate(main = ifelse(into == max(into), T, F))

    if (nrow(spk_pks_comp) > 1) {

      co_pid <- pull(spk_pks_comp, pid)
      spiked_pid <-  filter(spk_pks_comp, spiked == T) %>% pull(pid)

      ## get cor values between all co-eluting peaks, pair-wise
      # spk_pks_comp_cor <- expand.grid(pid1 = co_pid, pid2 = co_pid) %>%
      #   filter(pid1 != pid2) %>%
      #   group_by(pid1, pid2) %>%
      #   mutate(cor = (filter(cor_mat,
      #                        (x == pid1| x == pid2) &
      #                          (y == pid1| y == pid2)) %>% pull(cor)) )

      ## get cor values for alll co-eluting peaks
      cor_all <- cor_mat %>%
        ungroup() %>%
        filter(x != y) %>%
        filter(x %in% co_pid | y %in% co_pid) %>%
        select(pid1 = x, pid2 = y, cor) %>%
        mutate(comp = ifelse(pid1 %in% co_pid & pid2 %in% co_pid, T, F))

      ## extract cor values for between all co-eluting peaks, pair-wise
      cor_co <- cor_all %>%
        filter(comp == T)

      ## get correlations from real peak alignment
      cor_co_real <- buildCOR(co_ind = co_pid, eic = eic, pearson = pearson)

      ## extract EIC values for every corelated peak
      spec <- spk_pks_comp %>%
        group_by(pid) %>%
        do(extractSPECTRUM(co = ., raw = raw)) %>%
        ungroup()

      ## plot network of component's peaks
      buildNETWORK(poi_co_cor = cor_co_real,
                  pkscomps = pks,
                  co_ind = co_pid,
                  thr = thr,
                  pks = pks,
                  p = spiked_pid, # p is the pid of the main peak
                  plot = TRUE,
                  out_dir = out_dir,
                  fname = fname,
                  return = F)

      ####---- save MZ differences within the component
      cor_co_real <- cor_co_real %>%
        select(pid1 = x, pid2 = y, cor) %>%
        group_by(pid1, pid2) %>%
        mutate(pid_pid = paste(min(c(pid1, pid2)), max(c(pid1, pid2)), sep = "_")) %>%
        group_by(pid_pid) %>%
        slice(1) %>%
        mutate(mz_pid1 = pks %>% filter(pid == pid1) %>% pull(mz),
               mz_pid2 = pks %>% filter(pid == pid2) %>% pull(mz)) %>%
        mutate(mz_dif = round(abs(mz_pid1 - mz_pid2), digits = 3)) %>%
        ungroup()

      write.csv(cor_co_real, file = paste0(out_dir, "/", fname, "_spikedMR-", spk[i,"formula"], "_component_cor.csv"), row.names = F, quote = F)


      ####---- plot EICs of main adduct and all co-eluting peaks ----
      g_labels <- cor_co_real %>%
        filter(pid1 == spiked_pid | pid2 == spiked_pid) %>%
        mutate(not_spiked_pid = ifelse(pid1 == spiked_pid, pid2, pid1)) %>%
        group_by(not_spiked_pid) %>%
        slice(1) %>%
        mutate(g_text = paste0("Correlation: ", round(cor, digits = 3), ". MZ diff: ", mz_dif)) %>%
        ungroup() %>%
        select(pid1, pid2, not_spiked_pid, g_text) %>%
        bind_rows(data.frame(pid1 = spiked_pid, pid2 = spiked_pid, not_spiked_pid = spiked_pid, g_text = paste("Spiked MR peak PID: ", spiked_pid), stringsAsFactors = F)) %>%
        arrange(not_spiked_pid)

      g <- ggplot(spec, aes(x = rtime, xend = rtime, y = 0, yend = intensity)) +
        ## exclusion of column 'co_pid' allows to retain data points in each facet
        geom_segment(data = filter(spec, spiked == T) %>% select(., -co_pid),
                     alpha = 0.6, colour = "black", na.rm = TRUE) +
        geom_point(data = filter(spec, spiked == T)  %>% select(., -co_pid), aes(y = intensity),
                   alpha = 0.6, colour = "black", na.rm = TRUE) +
        geom_segment(data = filter(spec, spiked == F),
                     aes(colour = as.factor(co_pid)),
                     alpha = 0.8, na.rm = TRUE) +
        geom_point(data = filter(spec, spiked == F),
                   aes(y = intensity,colour = as.factor(co_pid)),
                   alpha = 0.8, na.rm = TRUE) +
        scale_colour_brewer(palette = "Paired",
                            name = "Peak PID") +
        facet_wrap(~co_pid,
                   labeller =  as_labeller(setNames(g_labels$g_text, nm = g_labels$not_spiked_pid)),
                   scales = "free") +
        geom_hline(aes(yintercept = paramCWT@prefilter[2],
                       linetype = paste0("Prefilter: c(", paramCWT@prefilter[1], "," , paramCWT@prefilter[2], ")")),
                   colour= "#1b7837") +
        geom_hline(aes(yintercept = paramCWT@noise,
                       linetype = paste("Noise:", paramCWT@noise)),
                   colour= '#762a83') +
          scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                                guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837")))) +
        xlab("Retention time") +
        ylab("Intensity") +
        theme_bw() +
        theme(legend.position = "bottom")

      ggsave(filename = paste0(out_dir, "/", fname, "_spikedMR-", spk[i,"formula"], "_component_spectra.png"), plot = g,width = 10, height = 8, units = "in")

      ####---- save MZ differences with the real compound molecular weight ----

      spk_pks_comp <- spk_pks_comp %>%
        mutate(mz_diff_mw = spk[i,"mol_weight"] - mz,
               mz_diff_main = pks %>% filter(pid == spiked_pid) %>% pull(mz) - mz)
      write.csv(spk_pks_comp, file = paste0(out_dir, "/", fname, "_spikedMR-", spk[i,"formula"], "_component_pks.csv"), row.names = F, quote = F)


      ####---- plot spiked MR component's peaks vs non-related peaks ----
      ## use artificially-aligned correlation values
      g <- ggplot(cor_all) +
        geom_histogram(aes(x = cor, group = comp, fill = comp),
                       na.rm = TRUE) +
        scale_fill_brewer(name = "", labels = c("Unrelated peaks", "Component's peak")) +
        facet_wrap(~comp,
                   scales = "free") +
        ylab("Count") +
        xlab("Correlation coefficient") +
        theme_bw() +
        theme(legend.position = "bottom")

      ggsave(filename = paste0(out_dir, "/", fname, "_spikedMR-", spk[i,"formula"], "_component_cor.png"), plot = g,width = 10, height = 8, units = "in")

    } else {print("No co-eluting peaks")}


  }


}

checkANYpeak <- function(p_list, fname, pks, pks_comps, raw, eic, out_dir, paramCWT, thr = 0.95) {


  for(p in p_list) {

    message("Checking peak with pno: ", p )

    ## get all co-eluting peaks for the peak-of-interest
    p_scpos <- c(pks[p, "scpos"] - 1, pks[p, "scpos"] + 1)
    p_co <- pks %>%
      filter(between(scpos, p_scpos[1], p_scpos[2])) %>%
      mutate(poi = ifelse(pno == p, T, F)) %>%
      ## which peak is the most intense - i.e. the main
      mutate(main = ifelse(into == max(into), T, F))

    if (nrow(p_co) > 1) {

      p_co_pno <- pull(p_co, pno)

      ## get cor values between all co-eluting peaks, pair-wise
      p_co_cor <- buildCOR(co_ind = p_co_pno, eic = eic, pearson = T)

      ## extract EIC values for every corelated peak
      spec <- p_co %>%
        group_by(pno) %>%
        do(extractSPECTRUM(co = ., raw = raw)) %>%
        ungroup()

      ##  save MZ differences within the component
      p_co_cor <- p_co_cor %>%
        group_by(x, y) %>%
        mutate(x_y = paste(min(c(x, y)), max(c(x, y)), sep = "_")) %>%
        group_by(x, y) %>%
        slice(1) %>%
        mutate(mz_x = pks %>% filter(pno == x) %>% pull(mz),
               mz_y = pks %>% filter(pno == y) %>% pull(mz)) %>%
        mutate(mz_dif = round(abs(mz_x - mz_y), digits = 2)) %>%
        ungroup()

      ####---- plot EICs of main adduct and all co-eluting peaks
      g_labels <- p_co_cor %>%
        filter(x == p | y == p) %>%
        mutate(not_p_pno = ifelse(x == p, y, x)) %>%
        group_by(not_p_pno) %>%
        slice(1) %>%
        mutate(g_text = paste0("Cor: ", round(cor, digits = 2), ". MZ diff: ", mz_dif)) %>%
        ungroup() %>%
        select(x, y, not_p_pno, g_text) %>%
        bind_rows(data.frame(x = p, y = p, not_p_pno = p, g_text = paste("Peak-of-Interest"), stringsAsFactors = F)) %>%
        arrange(not_p_pno)

      g <- ggplot(spec, aes(x = rtime, xend = rtime, y = 0, yend = intensity)) +
        ## exclusion of column 'co_pid' allows to retain data points in each facet
        geom_segment(data = filter(spec, main == T) %>% select(., -co_pno),
                     alpha = 0.6, colour = "black", na.rm = TRUE) +
        geom_point(data = filter(spec, main == T)  %>% select(., -co_pno), aes(y = intensity),
                   alpha = 0.6, colour = "black", na.rm = TRUE) +
        geom_segment(data = filter(spec, main == F),
                     aes(colour = as.factor(co_pno)),
                     alpha = 0.8, na.rm = TRUE) +
        geom_point(data = filter(spec, main == F),
                   aes(y = intensity,colour = as.factor(co_pno)),
                   alpha = 0.8, na.rm = TRUE) +
        scale_color_viridis_d( name = "Peak ID") +
        facet_wrap(~co_pno,
                   labeller =  as_labeller(setNames(g_labels$g_text, nm = g_labels$not_p_pno)),
                   scales = "free_y") +
        # geom_hline(aes(yintercept = paramCWT@prefilter[2],
        #                linetype = paste0("Prefilter: c(", paramCWT@prefilter[1], "," , paramCWT@prefilter[2], ")")),
        #            colour= "#1b7837") +
        # geom_hline(aes(yintercept = paramCWT@noise,
        #                linetype = paste("Noise:", paramCWT@noise)),
        #            colour= '#762a83') +
        # scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
        #                       guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837")))) +
        xlab("Retention time") +
        ylab("Intensity") +
        theme_bw() +
        theme(legend.position = "bottom")

        grid::grid.newpage()
        pdf(width = 12, height = 8,  paper = "a4",
            file = paste0(out_dir, "/", fname, "_poi-", p, "_allCO_spectra.pdf"))
        grid::grid.draw(ggplotGrob(g))
        dev.off()

        ####---- plot network of all co-eluting peaks
        buildNETWORK(poi_co_cor = p_co_cor,
                     pkscomps = pks,
                     co_ind = p_co_pno,
                     thr = thr,
                     pks = pks,
                     p = p, # p is the pid of the main peak
                     plot = TRUE,
                     out_dir = out_dir,
                     fname = fname,
                     return = F)

    } else {print("No co-eluting peaks")}

}
}










#' Title
#'
#' @param co
#' @param raw
#'
#' @return
#' @export
#'
#' @examples
extractSPECTRUM <- function(co, raw) {

  eic_co <- xcms::chromatogram(raw, aggregationFun = "sum",
                               rt = matrix(c(co$rtmin, co$rtmax), ncol = 2),
                               mz = matrix(c(co$mzmin, co$mzmax), ncol = 2))

  spec <- data.frame(rtime = as.numeric(xcms::rtime(eic_co[1])),
             intensity = as.numeric(xcms::intensity(eic_co[1]))) %>%
    mutate(co_pno = co$pno) %>%
    mutate(poi = co$poi) %>%
    mutate(main = co$main)

  return(spec)

}

