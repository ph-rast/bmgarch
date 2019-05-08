.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
  op = options()
  options(mc.cores = parallel::detectCores())
}
  
