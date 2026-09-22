# Audit F-BZ7-3: the Rachlin brms model is parameterised on a normalised delay
# (x / m, m = median positive delay) so the default logk prior -- and hence the
# curve prior -- does not depend on the delay unit. Rachlin mu = 1/(1 + k x^s);
# a unit change x_c = x / c gives k_c = c^s k, so a data-unit logk prior centred
# at -log(m) is unit-invariant only when s = 1. The sampler sees
# kn = k m^s; every reported draw is back-transformed to data units:
# log k = log kn - s log m (only the intercept shifts: s is population-level,
# factor contrasts and subject offsets are differences).

skip_if_not_installed("brms")
skip_on_ci()

rach_data <- function(scale = 1) {
  set.seed(5)
  d <- expand.grid(id = factor(1:6), x = c(0, 1, 7, 30, 90, 365, 730))
  d$y <- pmin(pmax(1 / (1 + 0.02 * d$x^0.8) + rnorm(nrow(d), 0, 0.03), 0.01), 0.99)
  d$x <- d$x / scale
  d
}

mock_draws <- function() {
  set.seed(1)
  n <- 50
  posterior::as_draws_matrix(data.frame(
    b_logk_Intercept = rnorm(n, 0.3, 0.2),
    b_logk_groupB = rnorm(n, -0.5, 0.1),
    b_logs_Intercept = rnorm(n, log(0.8), 0.1),
    r_id__logk.1.Intercept. = rnorm(n, 0, 0.3)
  ))
}

test_that("Rachlin default logk prior is identical for days and weeks (normalised delay)", {
  pd <- default_dd_priors("rachlin", family = "beta", data = rach_data(1))
  pw <- default_dd_priors("rachlin", family = "beta", data = rach_data(7))
  key <- pd$class == "b" & pd$coef == "Intercept" & pd$nlpar == "logk"
  expect_identical(pd$prior[key], "normal(0, 2.5)")
  strip <- function(p) { p <- as.data.frame(p); attr(p, "autoscale_info") <- NULL; p }
  expect_identical(strip(pd), strip(pw))
  m_d <- attr(pd, "autoscale_info")$delay_scale
  m_w <- attr(pw, "autoscale_info")$delay_scale
  expect_equal(m_d, stats::median(rach_data(1)$x[rach_data(1)$x > 0]))
  expect_equal(m_d / m_w, 7)
  # the other equations keep the data-unit anchor (they are unit-invariant)
  pg <- default_dd_priors("green-myerson", family = "beta", data = rach_data(1))
  keyg <- pg$class == "b" & pg$coef == "Intercept" & pg$nlpar == "logk"
  expect_identical(pg$prior[keyg],
                   paste0("normal(", beezdiscounting:::.dd_brms_fmt_num(-log(m_d)), ", 2.5)"))
  expect_null(attr(pg, "autoscale_info")$delay_scale)
})

test_that("Rachlin data columns are identical for days and weeks", {
  d1 <- rach_data(1); d7 <- rach_data(7)
  m1 <- stats::median(d1$x[d1$x > 0]); m7 <- stats::median(d7$x[d7$x > 0])
  c1 <- beezdiscounting:::.dd_brms_rachlin_cols(d1$x, m1)
  c7 <- beezdiscounting:::.dd_brms_rachlin_cols(d7$x, m7)
  expect_equal(c1, c7, tolerance = 1e-12)
  expect_identical(c1$xzero, as.numeric(d1$x == 0))
  expect_equal(c1$xsafe, ifelse(d1$x == 0, 1, d1$x / m1))
})

test_that("back-transform shifts only the logk intercept by -s log m", {
  dm <- mock_draws()
  out <- beezdiscounting:::.dd_brms_backtransform_logk(dm, 30)
  s <- exp(dm[, "b_logs_Intercept"])
  expect_equal(unclass(out[, "b_logk_Intercept"]),
               unclass(dm[, "b_logk_Intercept"] - s * log(30)))
  for (v in c("b_logk_groupB", "b_logs_Intercept", "r_id__logk.1.Intercept.")) {
    expect_identical(unclass(out[, v]), unclass(dm[, v]))
  }
  expect_identical(beezdiscounting:::.dd_brms_backtransform_logk(dm, 1), dm)
  expect_error(
    beezdiscounting:::.dd_brms_backtransform_logk(
      dm[, c("b_logk_groupB", "b_logs_Intercept")], 30),
    "intercept"
  )
})

test_that("the draws accessor back-transforms Rachlin fits; legacy fits are untouched", {
  dm <- mock_draws()
  obj <- list(brmsfit = dm, param_info = list(equation = "rachlin", delay_scale = 30))
  out <- beezdiscounting:::.dd_brms_draws_matrix(obj)
  expect_equal(unclass(out[, "b_logk_Intercept"]),
               unclass(dm[, "b_logk_Intercept"] -
                         exp(dm[, "b_logs_Intercept"]) * log(30)))
  # a fit saved before delay_scale existed (data-unit draws) is not shifted,
  # even though its autoscale_info carries a median delay
  legacy <- list(brmsfit = dm, param_info = list(equation = "rachlin"),
                 autoscale_info = list(median_delay = 30))
  expect_identical(beezdiscounting:::.dd_brms_draws_matrix(legacy), dm)
  gm <- list(brmsfit = dm, param_info = list(equation = "green-myerson",
                                             delay_scale = 30))
  expect_identical(beezdiscounting:::.dd_brms_draws_matrix(gm), dm)
})

test_that("shared normalised draws give unit-consistent k and identical curves", {
  dm <- mock_draws()
  m_day <- 90; m_week <- 90 / 7
  lk_day <- beezdiscounting:::.dd_brms_backtransform_logk(dm, m_day)[, "b_logk_Intercept"]
  lk_week <- beezdiscounting:::.dd_brms_backtransform_logk(dm, m_week)[, "b_logk_Intercept"]
  s <- exp(dm[, "b_logs_Intercept"])
  expect_equal(unclass(lk_week), unclass(lk_day + s * log(7)))
  x_day <- c(1, 7, 30, 365)
  mu_day <- sapply(x_day, function(x) 1 / (1 + exp(lk_day) * x^s))
  mu_week <- sapply(x_day / 7, function(x) 1 / (1 + exp(lk_week) * x^s))
  expect_equal(mu_day, mu_week, tolerance = 1e-12)
})

test_that("predict() newdata columns are always rebuilt from x on the fitted scale", {
  obj <- list(param_info = list(equation = "rachlin", delay_scale = 30))
  nd <- data.frame(x = c(0, 15, 60), xzero = 0, xsafe = c(0, 15, 60))
  out <- beezdiscounting:::.dd_brms_prep_newdata(obj, nd)
  expect_identical(out$xzero, c(1, 0, 0))
  expect_equal(out$xsafe, c(1, 0.5, 2))
  legacy <- list(param_info = list(equation = "rachlin"))
  expect_equal(beezdiscounting:::.dd_brms_prep_newdata(legacy, nd)$xsafe,
               c(1, 15, 60))
  gm <- list(param_info = list(equation = "mazur"))
  expect_identical(beezdiscounting:::.dd_brms_prep_newdata(gm, nd), nd)
})

test_that("diagnostics draws gain the data-unit logk intercept for Rachlin", {
  arr <- posterior::as_draws_array(mock_draws())
  out <- beezdiscounting:::.dd_brms_diag_draws(arr, 30)
  expect_true("b_logk_Intercept_data_units" %in% posterior::variables(out))
  expect_true(all(posterior::variables(arr) %in% posterior::variables(out)))
  expect_identical(beezdiscounting:::.dd_brms_diag_draws(arr, 1), arr)
})

test_that("days vs weeks give the same Rachlin curve posterior (sampling)", {
  skip_if_not(identical(Sys.getenv("BEEZ_RUN_BRMS_TESTS"), "true"),
              "Set BEEZ_RUN_BRMS_TESTS=true to run brms sampling tests")
  skip_on_cran()
  fit_u <- function(scale) {
    suppressWarnings(fit_dd_brms(
      rach_data(scale), equation = "rachlin", boundary = "squeeze",
      chains = 2, iter = 1500, seed = 11, init = "prior_center",
      loo = FALSE, verbose = 0
    ))
  }
  fd <- fit_u(1); fw <- fit_u(7)
  expect_equal(fd$param_info$delay_scale / fw$param_info$delay_scale, 7)
  # Monte-Carlo-aware comparison (separate sampling runs): log k shifts by
  # s log 7 and the population curve at matched delays agrees.
  s_d <- exp(fd$model$coefficients[["log_s"]])
  expect_equal(fw$model$coefficients[["beta_k"]],
               fd$model$coefficients[["beta_k"]] + s_d * log(7), tolerance = 0.05)
  nd_d <- data.frame(x = c(1, 7, 30, 365))
  nd_w <- data.frame(x = nd_d$x / 7)
  pd <- predict(fd, newdata = nd_d, level = "population")$.fitted
  pw <- predict(fw, newdata = nd_w, level = "population")$.fitted
  expect_equal(pd, pw, tolerance = 0.02)
})
