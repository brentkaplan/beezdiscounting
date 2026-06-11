# TICKET-039: brms tier default priors -------------------------------------------
#
# Verifies default_dd_priors() / default_dd_choice_priors() against the design
# (beezdemand internal_docs/design/DESIGN-brms-tier.md section 3.2).

skip_if_not_installed("brms")

fmt6 <- function(x) format(x, digits = 6, scientific = FALSE, trim = TRUE)

prior_row <- function(pri, class = "", coef = "", nlpar = "", dpar = "", group = "") {
  out <- pri[
    pri$class == class & pri$coef == coef & pri$nlpar == nlpar &
      pri$dpar == dpar & pri$group == group,
  ]
  expect_identical(nrow(out), 1L)
  out
}

dd_priors_data <- function() {
  data.frame(
    id = rep(1:4, each = 5),
    x = rep(c(1, 7, 30, 90, 365), times = 4),
    y = rep(c(0.95, 0.8, 0.6, 0.4, 0.2), times = 4)
  )
}

test_that("autoscaled logk centers k * median(delay) = 1", {
  d <- dd_priors_data()
  pri <- default_dd_priors("mazur", family = "beta", data = d)

  info <- attr(pri, "autoscale_info")
  expect_identical(info$median_delay, 30)

  expect_identical(
    prior_row(pri, "b", coef = "Intercept", nlpar = "logk")$prior,
    paste0("normal(", fmt6(-log(30)), ", 2.5)")
  )
  expect_identical(prior_row(pri, "sd", nlpar = "logk")$prior, "student_t(3, 0, 1)")
  expect_identical(prior_row(pri, "phi")$prior, "gamma(2, 0.1)")
})

test_that("static fallbacks and family-specific dispersion priors", {
  pri_beta <- default_dd_priors("mazur", family = "beta")
  expect_null(attr(pri_beta, "autoscale_info"))
  expect_identical(
    prior_row(pri_beta, "b", coef = "Intercept", nlpar = "logk")$prior,
    "normal(-4.5, 2.5)"
  )

  pri_gauss <- default_dd_priors("mazur", family = "gaussian")
  expect_identical(prior_row(pri_gauss, "sigma")$prior, "student_t(3, 0, 0.25)")
  expect_identical(nrow(pri_gauss[pri_gauss$class == "phi", ]), 0L)
})

test_that("logs prior exists only for two-parameter equations", {
  pri_gm <- default_dd_priors("green-myerson", family = "beta")
  expect_identical(
    prior_row(pri_gm, "b", coef = "Intercept", nlpar = "logs")$prior,
    "normal(0, 0.5)"
  )
  pri_mz <- default_dd_priors("mazur", family = "beta")
  expect_identical(nrow(pri_mz[pri_mz$nlpar == "logs", ]), 0L)
})

test_that("choice priors carry loggamma and optional b0", {
  pri <- default_dd_choice_priors("mazur")
  expect_identical(
    prior_row(pri, "b", coef = "Intercept", nlpar = "loggamma")$prior,
    "normal(1, 1)"
  )
  expect_identical(nrow(pri[pri$nlpar == "b0", ]), 0L)

  pri_int <- default_dd_choice_priors("mazur", intercept = TRUE)
  expect_identical(
    prior_row(pri_int, "b", coef = "Intercept", nlpar = "b0")$prior,
    "normal(0, 1.5)"
  )
})

test_that("factor designs add the fold-change coefficient prior in the ACCESSOR (Codex 039-R2)", {
  pri <- default_dd_priors("mazur", family = "beta", factors = "group")
  row <- pri[pri$class == "b" & pri$coef == "" & pri$nlpar == "logk", ]
  expect_identical(nrow(row), 1L)
  expect_identical(row$prior, "normal(0, 1)")

  pri0 <- default_dd_priors("mazur", family = "beta")
  expect_identical(
    nrow(pri0[pri0$class == "b" & pri0$coef == "" & pri0$nlpar == "logk", ]), 0L
  )
})

test_that("user priors override defaults on the full key", {
  defaults <- default_dd_priors("mazur", family = "beta")
  user <- brms::set_prior("normal(-3, 1)", class = "b", coef = "Intercept", nlpar = "logk")
  merged <- beezdiscounting:::.dd_brms_merge_priors(user, defaults)
  expect_identical(
    prior_row(merged, "b", coef = "Intercept", nlpar = "logk")$prior,
    "normal(-3, 1)"
  )
  key <- paste(merged$class, merged$coef, merged$group, merged$dpar, merged$nlpar)
  expect_false(any(duplicated(key)))
})

test_that("default priors validate against every equation x family model", {
  d <- dd_priors_data()
  for (eq in c("mazur", "exponential", "green-myerson", "rachlin")) {
    for (fam in c("beta", "gaussian")) {
      spec <- beezdiscounting:::.dd_brms_formula(equation = eq, family = fam)
      dd <- d
      if (eq == "rachlin") {
        dd$xzero <- as.numeric(dd$x == 0)
        dd$xsafe <- ifelse(dd$x == 0, 1, dd$x)
      }
      pri <- default_dd_priors(eq, family = fam, data = d)
      expect_no_warning(expect_no_error(
        brms::validate_prior(pri, formula = spec$formula, data = dd, family = spec$family)
      ))
    }
  }
})

test_that("choice priors validate against the choice model", {
  d <- data.frame(
    id = rep(1:4, each = 4), delay = rep(c(1, 7, 30, 90), 4),
    rel = 2, choice = rbinom(16, 1, 0.5)
  )
  spec <- beezdiscounting:::.dd_brms_choice_formula(equation = "mazur")
  pri <- default_dd_choice_priors("mazur")
  expect_no_warning(expect_no_error(
    brms::validate_prior(pri, formula = spec$formula, data = d, family = spec$family)
  ))
})

test_that("choice factor designs add the coefficient prior in the ACCESSOR (TICKET-048)", {
  pri <- default_dd_choice_priors("mazur", factors = "group")
  row <- pri[pri$class == "b" & pri$coef == "" & pri$nlpar == "logk", ]
  expect_identical(nrow(row), 1L)
  expect_identical(row$prior, "normal(0, 1)")

  # intercept-only design: no class-level coefficient row (brms warns-as-unused)
  pri0 <- default_dd_choice_priors("mazur")
  expect_identical(
    nrow(pri0[pri0$class == "b" & pri0$coef == "" & pri0$nlpar == "logk", ]), 0L
  )
})

test_that("choice factor priors validate against the factor formula (TICKET-048)", {
  set.seed(7)
  d <- data.frame(
    id = rep(1:6, each = 4), delay = rep(c(1, 7, 30, 90), 6),
    rel = 2, choice = rbinom(24, 1, 0.5),
    group = factor(rep(c("ctrl", "treat"), each = 12))
  )
  spec <- beezdiscounting:::.dd_brms_choice_formula(
    equation = "mazur", factors = "group", data = d
  )
  pri <- default_dd_choice_priors("mazur", data = d, factors = "group")
  expect_no_warning(expect_no_error(
    brms::validate_prior(pri, formula = spec$formula, data = d, family = spec$family)
  ))
})
