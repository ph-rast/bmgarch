library(bmgarch)
library(testthat)
options(mc.cores = 2)

data_small <- stocks[1:100, c("toyota", "nissan")]

# ── rstan backend – MCMC ───────────────────────────────────────────────────────

fit_rstan <- suppressWarnings(
  bmgarch(data             = data_small,
          parameterization = "CCC",
          standardize_data = TRUE,
          iterations       = 50,
          chains           = 2,
          backend          = "rstan",
          sampling_algorithm = "MCMC")
)

fit_rstan2 <- suppressWarnings(
  bmgarch(data             = data_small,
          parameterization = "CCC",
          P = 2, Q = 2,
          standardize_data = TRUE,
          iterations       = 50,
          chains           = 2,
          backend          = "rstan",
          sampling_algorithm = "MCMC")
)

test_that("rstan MCMC CCC returns bmgarch object", {
  expect_is(fit_rstan, "bmgarch")
  expect_equal(fit_rstan$sampling_algorithm, "MCMC")
  expect_equal(fit_rstan$backend, "rstan")
})

test_that("rstan MCMC CCC summary works", {
  expect_true(!is.null(summary(fit_rstan)$model_summary))
})

test_that("rstan MCMC CCC forecast works", {
  fc <- forecast(fit_rstan, ahead = 3)
  expect_is(fc, "forecast.bmgarch")
})

test_that("rstan MCMC CCC model_weights works and sums to 1", {
  mlist <- bmgarch_list(fit_rstan, fit_rstan2)
  mw    <- suppressWarnings(model_weights(mlist, L = 98))
  expect_is(mw, "model_weights")
  expect_equal(sum(mw$wts), 1, tolerance = 1e-6)
})

# ── cmdstanr backend – MCMC ───────────────────────────────────────────────────

skip_if_not_installed("cmdstanr")

fit_cmd_mcmc <- suppressWarnings(
  bmgarch(data             = data_small,
          parameterization = "CCC",
          standardize_data = TRUE,
          iterations       = 50,
          chains           = 2,
          backend          = "cmdstanr",
          sampling_algorithm = "MCMC",
          seed             = 123,
          refresh          = 0)
)

fit_cmd_mcmc2 <- suppressWarnings(
  bmgarch(data             = data_small,
          parameterization = "CCC",
          P = 2, Q = 2,
          standardize_data = TRUE,
          iterations       = 50,
          chains           = 2,
          backend          = "cmdstanr",
          sampling_algorithm = "MCMC",
          seed             = 123,
          refresh          = 0)
)

test_that("cmdstanr MCMC CCC returns bmgarch object", {
  expect_is(fit_cmd_mcmc, "bmgarch")
  expect_equal(fit_cmd_mcmc$sampling_algorithm, "MCMC")
  expect_equal(fit_cmd_mcmc$backend, "cmdstanr")
})

test_that("cmdstanr MCMC CCC summary works", {
  expect_true(!is.null(summary(fit_cmd_mcmc)$model_summary))
})

test_that("cmdstanr MCMC CCC forecast works", {
  fc <- forecast(fit_cmd_mcmc, ahead = 3)
  expect_is(fc, "forecast.bmgarch")
})

test_that("cmdstanr MCMC CCC model_weights works and sums to 1", {
  mlist <- bmgarch_list(fit_cmd_mcmc, fit_cmd_mcmc2)
  mw    <- suppressWarnings(model_weights(mlist, L = 98))
  expect_is(mw, "model_weights")
  expect_equal(sum(mw$wts), 1, tolerance = 1e-6)
})

# ── cmdstanr backend – VB ─────────────────────────────────────────────────────

fit_cmd_vb <- suppressWarnings(
  bmgarch(data             = data_small,
          parameterization = "CCC",
          standardize_data = TRUE,
          iterations       = 500,
          backend          = "cmdstanr",
          sampling_algorithm = "VB",
          seed             = 123,
          refresh          = 0)
)

fit_cmd_vb2 <- suppressWarnings(
  bmgarch(data             = data_small,
          parameterization = "CCC",
          P = 2, Q = 2,
          standardize_data = TRUE,
          iterations       = 500,
          backend          = "cmdstanr",
          sampling_algorithm = "VB",
          seed             = 123,
          refresh          = 0)
)

test_that("cmdstanr VB CCC returns bmgarch object", {
  expect_is(fit_cmd_vb, "bmgarch")
  expect_equal(fit_cmd_vb$sampling_algorithm, "VB")
  expect_equal(fit_cmd_vb$backend, "cmdstanr")
})

test_that("cmdstanr VB CCC summary works", {
  suppressWarnings(
    expect_true(!is.null(summary(fit_cmd_vb)$model_summary))
  )
})

test_that("cmdstanr VB CCC forecast works", {
  fc <- suppressWarnings(forecast(fit_cmd_vb, ahead = 3))
  expect_is(fc, "forecast.bmgarch")
})

test_that("cmdstanr VB CCC model_weights works and sums to 1", {
  mlist <- bmgarch_list(fit_cmd_vb, fit_cmd_vb2)
  mw    <- suppressWarnings(model_weights(mlist, L = 98))
  expect_is(mw, "model_weights")
  expect_equal(sum(mw$wts), 1, tolerance = 1e-6)
})
