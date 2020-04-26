# The setting of R_TESTS exists to work around an R bug.
# See https://github.com/hadley/testthat/issues/144
# This should be removed once the issue is resolved.
Sys.setenv("R_TESTS" = "")

library(testthat)
library(massFlowR)
library(xcms) # required to fix R CMD CHECK with xcms internal functions

test_check("massFlowR")
