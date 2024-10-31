.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)

  ## set core options
  op = options()
  options(mc.cores = parallel::detectCores())
  
  ## Check for optional presence of cmdstanr
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    packageStartupMessage(
      "The `cmdstanr` package is not installed. Some functionality may be limited. ",
      "To fully use `bmgarch`, consider installing `cmdstanr` from https://mc-stan.org/cmdstanr"
    )
  }
}
  
