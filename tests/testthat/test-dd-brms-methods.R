# TICKET-039: beezdiscounting_brms S3 methods --------------------------------
#
# Fixture-based (data-raw/make-dd-brms-fixtures.R). Mirrors the
# beezdiscounting_tmb contracts: the EXACT 8-column tidy table, dd glance
# columns, confint shape.

skip_if_not_installed("brms")
skip_on_ci()  # brms fits real Stan models; too slow/fragile under covr on CI (run locally)
skip_if_not_installed("posterior")

read_dd_fixture <- function() {
  path <- testthat::test_path("fixtures", "brms", "fit-mazur-beta.rds")
  skip_if_not(file.exists(path), "fixture missing: fit-mazur-beta")
  meta_path <- testthat::test_path("fixtures", "brms", "fixture-meta.rds")
  if (file.exists(meta_path)) {
    meta <- readRDS(meta_path)
    skip_if(
      utils::packageVersion("brms") < meta$brms_version,
      "installed brms older than fixture"
    )
  }
  readRDS(path)
}

test_that("tidy keeps the exact 8-column dd contract", {
  fit <- read_dd_fixture()
  t <- tidy(fit)

  expect_identical(
    names(t),
    c("term", "estimate", "std.error", "statistic", "p.value",
      "component", "estimate_scale", "term_display")
  )
  expect_true(all(is.na(t$statistic[t$component == "fixed"])))

  # behavioral: k natural estimate = median of exp draws
  draws <- posterior::as_draws_matrix(fit$brmsfit)
  kd <- as.numeric(draws[, "b_logk_Intercept"])
  k_row <- t[t$term == "k:(Intercept)", ]
  expect_equal(k_row$estimate, stats::median(exp(kd)))
  expect_equal(k_row$std.error, stats::sd(exp(kd)))
  expect_identical(k_row$estimate_scale, "natural")

  # variance rows follow the TMB reporting convention
  expect_true("sigma_u (log10-k RE SD)" %in% t$term)
  expect_true("phi (precision)" %in% t$term)
  expect_equal(
    t$estimate[t$term == "sigma_u (log10-k RE SD)"],
    stats::median(as.numeric(draws[, "sd_id__logk_Intercept"])) / log(10)
  )
})

test_that("tidy report_space internal/log10 transform the draws", {
  fit <- read_dd_fixture()
  draws <- posterior::as_draws_matrix(fit$brmsfit)
  kd <- as.numeric(draws[, "b_logk_Intercept"])

  ti <- tidy(fit, report_space = "internal", effects = "fixed")
  expect_equal(ti$estimate[ti$term == "k:(Intercept)"], stats::median(kd))

  t10 <- tidy(fit, report_space = "log10", effects = "fixed")
  expect_equal(
    t10$estimate[t10$term == "k:(Intercept)"],
    stats::median(kd / log(10))
  )
})

test_that("glance meets the dd contract with loo and diagnostics", {
  fit <- read_dd_fixture()
  g <- glance(fit)
  expect_identical(nrow(g), 1L)
  expect_true(all(c(
    "model_class", "backend", "equation", "family", "nobs", "n_subjects",
    "n_random_effects", "converged", "logLik", "AIC", "BIC",
    "elpd_loo", "looic", "rhat_max", "num_divergent"
  ) %in% names(g)))
  expect_identical(g$model_class, "beezdiscounting_brms")
  expect_identical(g$backend, "brms")
  expect_identical(g$logLik, NA_real_)
  expect_true(is.finite(g$elpd_loo))
})

test_that("confint returns quantile CIs on transformed draws", {
  fit <- read_dd_fixture()
  ci <- confint(fit)
  expect_true(all(c("term", "estimate", "conf.low", "conf.high", "level") %in% names(ci)))

  draws <- posterior::as_draws_matrix(fit$brmsfit)
  kn <- exp(as.numeric(draws[, "b_logk_Intercept"]))
  row <- ci[ci$term == "k:(Intercept)", ]
  expect_equal(row$conf.low, unname(stats::quantile(kn, 0.025)))
  expect_equal(row$conf.high, unname(stats::quantile(kn, 0.975)))

  expect_identical(nrow(confint(fit, parm = "beta_k")), 1L)
})

test_that("predict/fitted/residuals/augment behave", {
  fit <- read_dd_fixture()
  pr <- predict(fit)
  expect_identical(nrow(pr), fit$param_info$n_obs)
  ep <- brms::posterior_epred(fit$brmsfit, re_formula = NULL)
  expect_equal(pr$.fitted, apply(ep, 2, stats::median))

  expect_identical(predict(fit, type = "parameters"), fit$subject_pars)

  fv <- fitted(fit)
  expect_equal(residuals(fit), fit$data$y - fv)
  au <- augment(fit)
  expect_true(all(c(".fitted", ".resid", ".std_resid") %in% names(au)))
})

test_that("ranef, coef, nobs, logLik guard", {
  fit <- read_dd_fixture()
  re <- ranef(fit)
  expect_true(all(c("id", "r_logk", "k") %in% names(re)))
  expect_identical(coef(fit), fit$model$coefficients)
  expect_identical(nobs(fit), fit$param_info$n_obs)
  expect_error(logLik(fit), "loo")
})

test_that("print and summary render", {
  fit <- read_dd_fixture()
  expect_output(print(fit), "Discounting")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.beezdiscounting_brms")
  expect_output(print(sm), "Bayesian")
})
