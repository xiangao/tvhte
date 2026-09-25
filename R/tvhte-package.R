#' @keywords internal
"_PACKAGE"

## Base R functions used in package code. R CMD check reports these as
## undefined globals when they are not imported, because a package's own code
## does not see the attached search path.
#' @importFrom stats coef lm optim residuals rnorm var
NULL
