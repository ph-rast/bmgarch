library(bmgarch)
library(testthat)
options(mc.cores = 2)

# Minimal fast fit used across multiple tests
fit_const <- suppressWarnings(
  bmgarch(data         = stocks[1:100, c("toyota", "nissan")],
          parameterization = "const",
          standardize_data = TRUE,
          iterations   = 10)
)

# ── Basic object ───────────────────────────────────────────────────────────────

test_that("const returns a bmgarch object", {
  expect_is(fit_const, "bmgarch")
})

test_that("const param stored correctly", {
  expect_equal(fit_const$param, "const")
})

# ── H is truly constant across time ───────────────────────────────────────────

test_that("H_const is identical at every time step", {
  H_array <- rstan::extract(fit_const$model_fit, pars = "corH")[[1]]
  # H_array dims: [draws, T, nt, nt]
  # All T slices should be equal within each draw
  H_t1 <- H_array[, 1, , ]
  H_t2 <- H_array[, 2, , ]
  H_tT <- H_array[, dim(H_array)[2], , ]
  expect_equal(H_t1, H_t2)
  expect_equal(H_t1, H_tT)
})

# ── Mean structures ────────────────────────────────────────────────────────────

test_that("const with constant mean structure fits", {
  fit <- suppressWarnings(
    bmgarch(data         = stocks[1:100, c("toyota", "nissan")],
            parameterization = "const",
            meanstructure = "constant",
            standardize_data = TRUE,
            iterations   = 10)
  )
  expect_is(fit, "bmgarch")
  expect_equal(fit$meanstructure, 0L)
})

test_that("const with ARMA mean structure fits", {
  fit <- suppressWarnings(
    bmgarch(data         = stocks[1:100, c("toyota", "nissan")],
            parameterization = "const",
            meanstructure = "ARMA",
            standardize_data = TRUE,
            iterations   = 10)
  )
  expect_is(fit, "bmgarch")
  expect_equal(fit$meanstructure, 1L)
})

test_that("const with VAR mean structure fits", {
  fit <- suppressWarnings(
    bmgarch(data         = stocks[1:100, c("toyota", "nissan")],
            parameterization = "const",
            meanstructure = "VAR",
            standardize_data = TRUE,
            iterations   = 10)
  )
  expect_is(fit, "bmgarch")
  expect_equal(fit$meanstructure, 2L)
})

# ── Higher GARCH order (P,Q unused but must not error) ────────────────────────

test_that("const with P=2, Q=2 fits without error", {
  fit <- suppressWarnings(
    bmgarch(data         = stocks[1:100, c("toyota", "nissan")],
            parameterization = "const",
            P = 2, Q = 2,
            standardize_data = TRUE,
            iterations   = 10)
  )
  expect_is(fit, "bmgarch")
})

# ── Forecasting ────────────────────────────────────────────────────────────────

test_that("forecast returns forecast.bmgarch for const", {
  fc <- forecast(fit_const, ahead = 3)
  expect_is(fc, "forecast.bmgarch")
})

test_that("forecast H is constant across horizon for const", {
  fc <- forecast(fit_const, ahead = 5)
  # Posterior-mean H at step 1 and step 5 should be identical
  H_step1 <- fc$forecast$var[1, "mean", ]
  H_step5 <- fc$forecast$var[5, "mean", ]
  expect_equal(H_step1, H_step5)
})

# ── LFO-CV ────────────────────────────────────────────────────────────────────

test_that("loo (backward) works for const", {
  lfob <- suppressWarnings(loo(fit_const, mode = "backward", L = 99))
  expect_is(lfob, "loo.bmgarch")
})

# ── bmgarch_list and model_weights (const as baseline) ────────────────────────

test_that("const can be combined with MGARCH models in bmgarch_list", {
  fit_ccc <- suppressWarnings(
    bmgarch(data         = stocks[1:100, c("toyota", "nissan")],
            parameterization = "CCC",
            standardize_data = TRUE,
            iterations   = 10)
  )
  mlist <- bmgarch_list(fit_const, fit_ccc)
  expect_is(mlist, "bmgarch_list")
})

test_that("model_weights works with const in the candidate set", {
  fit_ccc <- suppressWarnings(
    bmgarch(data         = stocks[1:100, c("toyota", "nissan")],
            parameterization = "CCC",
            standardize_data = TRUE,
            iterations   = 10)
  )
  mlist <- bmgarch_list(fit_const, fit_ccc)
  mw    <- suppressWarnings(model_weights(mlist, L = 98))
  expect_is(mw, "model_weights")
  # Weights must sum to 1
  expect_equal(sum(mw$wts), 1, tolerance = 1e-6)
})

# ── Summary / print / plot ─────────────────────────────────────────────────────

test_that("summary returns non-null model_summary for const", {
  expect_true(!is.null(summary(fit_const)$model_summary))
})

test_that("print returns bmgarch object for const", {
  expect_is(print(fit_const)$RTS_full, "array")
})

test_that("plot returns list of ggplots for const", {
  out <- plot(fit_const, askNewPage = FALSE)
  expect_is(out, "list")
  dev.off( )
})
