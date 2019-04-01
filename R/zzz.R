.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste("\nThis is massFlowR version", packageVersion("massFlowR"), ": development", "\n"))
}

.onLoad <- function(libname = find.package("massFlowR"), pkgname = "massFlowR"){
  ## quiets R CMD check notes about the variables used inside foreach loops, and upon loading RDA files (chem.file)
  if(getRversion() >= "2.15.1") {
    utils::globalVariables(c("f", "p", "s", "n", "pkg", "sn", "peakid", "chem.file")) 
    invisible()
  }
}
