.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0(
      "\nThis is massFlowR version ", utils::packageVersion("massFlowR"), ": development", "\n",
      "Date of built: ", utils::packageDate("massFlowR")
    )
  )
}

.onLoad <- function(libname = find.package("massFlowR"), pkgname = "massFlowR") {
  ## quiets R CMD check notes about the variables used inside foreach loops, and
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
      ## variables used by functions with foreach backend
      "f", "p", "s", "n", "pkg", "sn", "peakid",
      ## variables used upon loading RDA files by buildDB
      "chem.file",
      # data.frame column names used in ggplot (cannot use aes_string due to transformations to the column in aes())
      "color_by", "into", "mz", "run_order"
    ))
    invisible()
  }
}
