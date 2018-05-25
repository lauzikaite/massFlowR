#' Title
#'
#' @param raw
#' @param pks
#'
#' @return
#' @export
#'
#' @examples
extractEIC <- function(raw = raw, pks = pks) {

  start <- Sys.time()
  eic <- xcms::chromatogram(raw,
                            rt = data.frame(
                              rt_lower = pks$rtmin,
                              rt_upper = pks$rtmax),
                            mz = data.frame(
                              mz_lower = pks$mzmin,
                              mz_upper = pks$mzmax))
  clean_eic <- lapply(1:nrow(eic), function(ch) {
    clean(eic[ch, ], na.rm = T)
  })
  etime <- Sys.time()
  print(Sys.time() - start)
  return(clean_eic)
}
