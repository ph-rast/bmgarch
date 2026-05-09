# bmgarch (development version)

## cmdstanr backend support

* `bmgarch()` now accepts `backend = "cmdstanr"` in addition to the default
  `"rstan"`. Both MCMC (NUTS) and variational Bayes (`sampling_algorithm = "VB"`)
  are supported with the new backend.

* `summary()`, `print()`, `plot()`, `forecast()`, `fitted()`, `loo()`, and
  `model_weights()` all work with cmdstanr-fitted objects.

* Stan function `a_b_scale_jacobian` renamed to `a_b_scale_log_jac` to avoid a
  future name collision with Stan 2.39's reserved `_jacobian` suffix (affects
  `inst/stan/functions/jacobian.stan`).

* The `iterations` argument now defaults to `2000` for MCMC and `30000` for VB
  when left unspecified. Previously the single default of `2000` was too low for
  ADVI gradient ascent.

* The `backend` field is stored in every `bmgarch` object so that `loo()`
  refits and other operations that call back into `bmgarch()` automatically use
  the same backend.

* Threading defaults: cmdstanr models compiled with `stan_threads = TRUE` now
  default to `threads_per_chain = 1` (MCMC) and `threads = 1` (VB). Users can
  override via `...`.

## Bug fixes

* Fixed ensemble covariance underestimation in `.weighted_samples()`. The
  between-model variance term required by the law of total variance,
  `Σ_m w_m (μ̂_m − μ̂_mix)(μ̂_m − μ̂_mix)ᵀ`, is now added to the
  weighted combination of `H` and `H_forecasted`. This corrects the
  ensemble posterior predictive covariance for both backcasts and forecasts
  whenever models with differing mean structures are combined. The fix is a
  no-op for single-model objects.

* `model_weights()` and `loo()` now work with cmdstanr objects: `log_lik`
  extraction and chain/iteration metadata no longer rely on rstan S4 slot
  access (`@sim`).

* Weighted forecasts (`forecast(..., weights = mw)`) now work with cmdstanr
  backends via a new `.extract_param_list()` helper that replicates
  `rstan::extract()` output for cmdstanr draw objects.

* Forecast draws passed to `rstan::gqs()` are now correctly extracted via
  `posterior::as_draws_matrix()` for cmdstanr fits, fixing a "no dimnames"
  error.

## Tests

* Added `tests/testthat/test-backends.R` covering `summary`, `forecast`, and
  `model_weights` for all three combinations: rstan/MCMC, cmdstanr/MCMC, and
  cmdstanr/VB on a CCC model.

# bmgarch 2.0.0

* All stan models are rewritten to match the new `rstan 2.26.0` array syntax.

# bmgarch 1.1.0

* Meanstructure now also takes the `VAR` argument for a VAR(1) model. 
* Added variational Bayes `VB` as a sampling algorithm option. bmgarch now takes `VB` and  the standard `MCMC` argument in sampling_algorithm. `VB` is inherited from `rstan` and is still experimental - use with caution.
* Updated rstantools to facilitate rstan upgrade (according to rstan developer guide)
* standat now checks for constant vectors in data, and returns error if there's no variance (cf. discourse [mc-stan](https://discourse.mc-stan.org/t/getting-rejected-initial-values-running-bmgarch-in-r/24002))

# bmgarch 1.0.1

* `loo` function now takes `m`-ahead argument for the "backward" approximation allowing one to tune the forcast to arbitrary steps ahead.

* Fixed package dependency versions.

# bmgarch 1.0.0
Initial CRAN release
