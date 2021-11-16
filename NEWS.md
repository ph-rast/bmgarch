# bmgarch 1.1.0

* Meanstructure now also takes the `VAR` argument for a VAR(1) model. 
* Added variational Bayes `VB` as a sampling algorithm option. bmgarch now takes `VB` and  the standard `MCMC` argument in sampling_algorithm. `VB` is inherited from `rstan` and is still experimental - use with caution.
* Updated rstantools to facilitate rstan upgrade (according to rastan developer guide)
* standat now checks for constant vectors in data, and returns error

# bmgarch 1.0.1

* `loo` function now takes `m`-ahead argument for the "backward" approximation allowing one to tune the forcast to arbitrary steps ahead.

* Fixed package dependency versions.

# bmgarch 1.0.0
Initial CRAN release
