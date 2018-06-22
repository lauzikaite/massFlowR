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
    cmp <- paste0(spk[i,"name"], ": peak was NOT picked")
    g <- ggplot(gdf) +
      geom_segment(data = gdf,  aes(x = rtime, xend = rtime, y = 0, yend = intensity), alpha = 0.8, colour = "#a6cee3", na.rm = TRUE) +
      geom_point(data = gdf, aes(x = rtime, y = intensity), colour = "#a6cee3", na.rm = TRUE) +
      geom_hline(aes(yintercept= prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
      geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
      scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                            guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837")))) +
      ylab("intensity") +
      xlab("retention time")

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
        geom_segment(data = gdf,  aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(peak)), alpha = 0.8, na.rm = TRUE) +
        geom_point(data = gdf, aes(x = rtime, y = intensity, colour = as.factor(peak)), alpha = 0.8, na.rm = TRUE) +
        scale_colour_brewer(palette="Paired",
                          name = "", labels = c("Full spectrum", "Picked peak")) +
        geom_hline(aes(yintercept= prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
        geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
        scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                              guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837")))) +
        ylab("intensity") +
        xlab("retention time")

      ## save output
      match <- data.frame(sample_id = fname,
                          spiked_id = spk[i,"id"],
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

        match <- data.frame(sample_id = fname,
                            spiked_id = spk[i,"id"],
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
          geom_segment(data = gdf,  aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(peak)), alpha = 0.8, na.rm = TRUE) +
          geom_point(data = gdf, aes(x = rtime, y = intensity, colour = as.factor(peak)), alpha = 0.8, na.rm = TRUE) +
          scale_colour_brewer(palette="Paired",
                              name = "", labels = c("Full spectrum", "Picked peak")) +
          geom_hline(aes(yintercept= prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
          geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
          scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                                guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837")))) +
          ylab("intensity") +
          xlab("retention time")

        match <- data.frame(sample_id = fname,
                            spiked_id = spk[i,"id"],
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
      ggtitle(cmp) +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(filename = paste0(out_dir, "/", fname, "_spiked_", spk[i,"id"], ".png"), plot = g,width = 10, height = 8, units = "in")

    match_df <- rbind(match_df, match)
  }
  return(match_df)
}


checkSPIKEDcomps <- function(i, fname, pks, raw, eic, spk_pks, cor_mat, out_dir, prefilter, noise, rt_err = 3, mz_err = 0.005) {

  if (missing(i)) {

    ## for all spiked MR compounds
    i_list <- 1:nrow(spk_pks)

  } else {

    ## for selected spiked MR compound only
    i_list <- i
  }

  match_df <- data.frame()

  for(i in i_list) {

    ## get all peaks in the same component as the matched peak
    spk_pks_comp <- pks %>%
      filter(comp == (pks %>% filter(pid == spk_pks[i,"pks_id"]) %>% select(comp) %>% pull) ) %>%
      mutate(spiked = ifelse(pid == spk_pks[i,"pks_id"], T, F)) %>%
      ## which peak is the most intense - i.e. the main
      mutate(main = ifelse(into == max(into), T, F))

    if (nrow(spk_pks_comp) > 1) {

      co_pid <- pull(spk_pks_comp, pid)
      spiked_pid <-  filter(spk_pks_comp, spiked == T) %>% pull(pid)

      ## extract EIC values for every corelated peak
      spec <- spk_pks_comp %>%
        group_by(pid) %>%
        do(extractSPECTRUM(co = ., raw = raw))

      ## get cor values between all co-eluting peaks, pair-wise

      ## doesn't work yet since previous cor mat generations were wrong
      # spk_pks_comp_cor <- expand.grid(pid1 = co_pid, pid2 = co_pid) %>%
      #   filter(pid1 != pid2) %>%
      #   group_by(pid1, pid2) %>%
      #   mutate(cor =  (filter(cor_mat,
      #                        (x == pid1| x == pid2) &
      #                          (y == pid1| y == pid2)) %>% pull(cor)) )


      # ## pids won't match order in the peak table any more if reduced pks-comps-cls.txt file is used
      # co_ind <- which(pks$pid %in% co_pid)
      # spk_pks_comp_cor <- buildCOR(co_ind = co_ind, eic = eic, pearson = pearson)

      ## make spiked - coeluting pairs
      comp_cor <-  spk_pks_comp_cor %>%
        filter(x == spiked_pid | y == spiked_pid)


      ####---- plot EICs of main adduct and all co-eluting peaks ----
      ## version (A) - all together
      # g <- ggplot(spec) +
      #   geom_segment(data = spec,  aes(x = rtime, xend = rtime, y = 0, yend = intensity, colour = as.factor(peak)), alpha = 0.8) +
      #   geom_point(data = spec, aes(x = rtime, y = intensity, colour = as.factor(peak)), alpha = 0.8) +
      #   geom_hline(aes(yintercept= prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
      #   geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
      #   scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
      #                         guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837")))) +
      #   scale_colour_brewer(palette="Paired",
      #                       name = "") +
      #   ylab("intensity") +
      #   xlab("retention time")

      ## version B - facet



      g_text <- paste0("Correlation: ", round(unique(spec$cor), digits = 3), ". MZ difference: ", abs(round(unique(gdf$mz_diff), digits = 3)))


      g <- ggplot(spec, aes(x = rtime, xend = rtime, y = 0, yend = intensity)) +
        ## exclusion of column 'co_pid' allows to retain data points in each facet
        geom_segment(data = filter(spec, spiked == T) %>% select(., -co_pid),
                     alpha = 0.6, colour = "black") +
        geom_point(data = filter(spec, spiked == T)  %>% select(., -co_pid), aes(y = intensity),
                   alpha = 0.6, colour = "black") +
        geom_segment(data = filter(spec, spiked == F),
                     alpha = 0.6, aes(colour = as.factor(co_pid))) +
        scale_colour_brewer(palette = "Paired",
                            name = "Peak PID") +
        geom_point(data = filter(spec, spiked == F),
                   aes(y = intensity,colour = as.factor(co_pid)), alpha = 0.6) +
        facet_wrap(~co_pid) +
        geom_hline(aes(yintercept = prefilter[2], linetype = paste0("Prefilter: c(", prefilter[1], "," , prefilter[2], ")")), colour= "#1b7837") +
          geom_hline(aes(yintercept= noise, linetype = paste("Noise:", noise)), colour= '#762a83') +
          scale_linetype_manual(name = "centWave parameters", values = c(2, 2),
                                guide = guide_legend(override.aes = list(color = c("#762a83", "#1b7837")))) +
        xlab("Retention time") +
        ylab("Intensity")

      ggsave(filename = paste0(out_dir_fname, "/", fname, "_spiked_", spk[i,"id"], "_coelutingPP.png"), plot = g,width = 10, height = 8, units = "in")


      ####---- plot mz diff vs cor ----
      cmp <-  paste0(spk[i,"name"])
      g <- ggplot(data = all_poi_co) +
        geom_point(aes(x = cor, y = mz_diff)) +
        geom_text(aes(x = cor, y = mz_diff, label = round(mz_diff, digits = 3)), vjust = 0, nudge_y = 0.5) +
        scale_x_continuous(breaks = c(seq(-1, 1, by = 0.1))) +
        ylab("MZ difference") +
        xlab("Correlation coefficient") +
        theme_bw() +
        ggtitle(cmp) +
        theme(plot.title = element_text(hjust = 0.5))
      ggsave(filename = paste0(out_dir_fname, "/", fname, "_spiked_", spk[i,"id"], "_coelutingPP-mzdiff.png"), plot = g,width = 10, height = 8, units = "in")
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
    mutate(co_pid = co$pid) %>%
    mutate(spiked = co$spiked) %>%
    mutate(main = co$main)

  return(spec)

}

