##' Internal function to obtain stan code paths
##' @title Path to stan files
##' @return Full path to folder containing stan files
##' @author philippe
##' @keywords internal
.get_target_stan_path <- function() {
  package_path <-
    system.file(package = "bmgarch", "stan")
  return(package_path)
}
