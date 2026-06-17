# TICKET-039: fit_dd_brms() smoke tests ------------------------------------------
#
# Real sampling is opt-in (BEEZ_RUN_BRMS_TESTS=true); validation tests run
# anywhere brms is installed (they error before compilation).

skip_if_not_installed("brms")
skip_on_ci()  # brms fits real Stan models; too slow/fragile under covr on CI (run locally)

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

choice_smoke_data <- function(n_id = 10) {
  set.seed(41)
  d <- expand.grid(
    id = factor(seq_len(n_id)),
    delay = c(1, 7, 30, 90, 180),
    rep = 1:4
  )
  d$ss_amount <- 50
  d$ll_amount <- 100
  k_i <- exp(log(0.02) + rnorm(n_id, 0, 0.4))
  D <- 1 / (1 + k_i[d$id] * d$delay)
  eta <- 3 * ((d$ll_amount / d$ss_amount) * D - 1)
  d$choice <- rbinom(nrow(d), 1, plogis(eta))
  d$rep <- NULL
  d
}

test_that("fit_dd_choice_brms rejects descriptive mode with guidance", {
  d <- choice_smoke_data()
  expect_error(
    fit_dd_choice_brms(d, mode = "descriptive"),
    "descriptive"
  )
})

test_that("fit_dd_choice_brms structural: object contract and recovery", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("BEEZ_RUN_BRMS_TESTS"), "true"),
    "Set BEEZ_RUN_BRMS_TESTS=true to run brms sampling tests"
  )

  d <- choice_smoke_data()
  fit <- fit_dd_choice_brms(
    d,
    equation = "mazur",
    chains = 2, iter = 600, warmup = 300,
    cores = 2, seed = 123,
    loo = FALSE, verbose = 0
  )

  expect_s3_class(fit, "beezdiscounting_choice_brms")
  expect_identical(names(fit$model$coefficients), c("beta_k", "log_gamma"))
  # truth: log k = log(0.02), gamma = 3
  expect_lt(abs(fit$model$coefficients[["beta_k"]] - log(0.02)), 1)
  expect_lt(abs(exp(fit$model$coefficients[["log_gamma"]]) - 3), 2.5)

  t <- tidy(fit)
  expect_identical(
    names(t),
    c("term", "estimate", "std.error", "statistic", "p.value",
      "component", "estimate_scale", "term_display")
  )
  expect_true(all(c("k:(Intercept)", "gamma") %in% t$term))
  # TMB choice contract: k is the fixed row; gamma is a shape row (Codex 039-B2)
  expect_identical(t$component[t$term == "k:(Intercept)"], "fixed")
  expect_identical(t$component[t$term == "gamma"], "shape")
  expect_identical(t$term_display[t$term == "gamma"], "gamma")

  # confint with parm aliases (Codex 039-B1)
  ci <- confint(fit)
  expect_true(all(c("term", "estimate", "conf.low", "conf.high", "level") %in% names(ci)))
  expect_identical(nrow(confint(fit, parm = "log_gamma")), 1L)
  expect_identical(nrow(confint(fit, parm = "k:(Intercept)")), 1L)

  pr <- predict(fit)
  expect_true(all(pr$.fitted >= 0 & pr$.fitted <= 1))
  expect_identical(nrow(fit$subject_pars), 10L)

  # intercept-only EMMs no longer abort (TICKET-048): population k row
  em <- get_dd_param_emms(fit)
  expect_identical(em$level, "(Intercept)")
  expect_equal(em$k_log, fit$model$coefficients[["beta_k"]])
})

test_that("fit_dd_choice_brms recovers a between-subject k effect (TICKET-048)", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("BEEZ_RUN_BRMS_TESTS"), "true"),
    "Set BEEZ_RUN_BRMS_TESTS=true to run brms sampling tests"
  )

  # two groups of 8 subjects; truth: log k_treat - log k_ctrl = log(4)
  set.seed(43)
  n_per_group <- 8
  d <- expand.grid(
    id = seq_len(2 * n_per_group),
    delay = c(1, 7, 30, 90, 180),
    rep = 1:4
  )
  d$group <- factor(ifelse(d$id <= n_per_group, "ctrl", "treat"))
  d$ss_amount <- 50
  d$ll_amount <- 100
  logk_i <- log(0.01) + log(4) * (seq_len(2 * n_per_group) > n_per_group) +
    rnorm(2 * n_per_group, 0, 0.3)
  k_i <- exp(logk_i)
  D <- 1 / (1 + k_i[d$id] * d$delay)
  eta <- 3 * ((d$ll_amount / d$ss_amount) * D - 1)
  d$choice <- rbinom(nrow(d), 1, plogis(eta))
  d$id <- factor(d$id)
  d$rep <- NULL

  fit <- fit_dd_choice_brms(
    d,
    equation = "mazur", factors = "group",
    chains = 2, iter = 600, warmup = 300,
    cores = 2, seed = 123,
    loo = FALSE, verbose = 0
  )

  coefs <- fit$model$coefficients
  expect_identical(sum(names(coefs) == "beta_k"), 2L)
  k_coefs <- unname(coefs[names(coefs) == "beta_k"])
  # acceptance: recovered group effect near log(4) ~ 1.386
  expect_lt(abs(k_coefs[2] - log(4)), 0.8)
  expect_lt(abs(k_coefs[1] - log(0.01)), 1)

  # EMM/contrast surface works end-to-end on a sampled factor fit
  em <- get_dd_param_emms(fit)
  expect_identical(nrow(em), 2L)
  cmp <- get_dd_comparisons(fit)
  expect_identical(nrow(cmp$k$contrasts_log10), 1L)
  expect_true(cmp$k$contrasts_ratio$ratio != 1)
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
