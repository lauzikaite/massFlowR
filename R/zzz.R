.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste("\nThis is massFlowR version", packageVersion("massFlowR"), "\n"))
}
