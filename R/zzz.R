.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
  op = options()
  options(mc.cores = parallel::detectCores())
    ## Check for presence of cmdstanr
  if( !require("cmdstanr") ) stop("The `cmdstanr` library as well as CmdStan are required to run this pacakge. Check out https://mc-stan.org/cmdstanr")
}
  
