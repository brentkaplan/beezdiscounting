# TICKET-039: fit_dd_brms() smoke tests ------------------------------------------
#
# Real sampling is opt-in (BEEZ_RUN_BRMS_TESTS=true); validation tests run
# anywhere brms is installed (they error before compilation).

skip_if_not_installed("brms")

dd_smoke_data <- function(n_id = 8) {
  set.seed(31)
  delays <- c(1, 7, 30, 90, 180, 365)
  d <- expand.grid(id = factor(seq_len(n_id)), x = delays)
  k_i <- exp(log(0.02) + rnorm(n_id, 0, 0.4))
  mu <- 1 / (1 + k_i[d$id] * d$x)
  d$y <- pmin(pmax(mu + rnorm(nrow(d), 0, 0.05), 0), 1)
  d
}

# --- validation (no sampling) ----------------------------------------------------

test_that("fit_dd_brms rejects k + phi random effects in v1", {
  d <- dd_smoke_data()
  expect_error(
    fit_dd_brms(d, random_effects = k + phi ~ 1),
    "k \\+ phi|phi random effects"
  )
})

test_that("fit_dd_brms sltb errors with the beta-analog pointer", {
  d <- dd_smoke_data()
  expect_error(fit_dd_brms(d, family = "sltb"), "beta")
})

test_that("boundary = 'error' refuses boundary responses", {
  d <- dd_smoke_data()
  d$y[1] <- 1
  expect_error(
    fit_dd_brms(d, family = "beta", boundary = "error"),
    "boundary"
  )
})

# --- sampling (opt-in) -------------------------------------------------------------

test_that("fit_dd_brms mazur/beta: object contract and recovery", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("BEEZ_RUN_BRMS_TESTS"), "true"),
    "Set BEEZ_RUN_BRMS_TESTS=true to run brms sampling tests"
  )

  d <- dd_smoke_data()
  fit <- fit_dd_brms(
    d,
    equation = "mazur", family = "beta",
    chains = 2, iter = 600, warmup = 300,
    cores = 2, seed = 123,
    loo = FALSE, verbose = 0
  )

  expect_s3_class(fit, "beezdiscounting_brms")
  expect_s3_class(fit$brmsfit, "brmsfit")
  expect_identical(names(fit$model$coefficients), "beta_k")
  # truth: log k = log(0.02) ~= -3.9
  expect_lt(abs(fit$model$coefficients[["beta_k"]] - log(0.02)), 1)

  vc <- fit$model$variance_components
  expect_true(all(c("sigma_u (log10-k RE SD)", "phi (precision)") %in% vc$Component))

  sp <- fit$subject_pars
  expect_identical(nrow(sp), 8L)
  expect_true(all(c("id", "k", "k_lower", "k_upper") %in% names(sp)))
  expect_true(all(sp$k > 0))

  expect_true(is.logical(fit$converged))
  expect_identical(fit$mcmc_info$backend, "rstan")
  expect_identical(fit$loglik, NA_real_)
  expect_identical(fit$AIC, NA_real_)
})

test_that("fit_dd_brms gaussian green-myerson fits with the s parameter", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("BEEZ_RUN_BRMS_TESTS"), "true"),
    "Set BEEZ_RUN_BRMS_TESTS=true to run brms sampling tests"
  )

  d <- dd_smoke_data()
  fit <- fit_dd_brms(
    d,
    equation = "green-myerson", family = "gaussian",
    chains = 2, iter = 600, warmup = 300,
    cores = 2, seed = 123,
    loo = TRUE, verbose = 0
  )

  expect_identical(names(fit$model$coefficients), c("beta_k", "log_s"))
  expect_s3_class(fit$loo, "loo")
  expect_true("sigma_e (Residual SD)" %in% fit$model$variance_components$Component)
})
