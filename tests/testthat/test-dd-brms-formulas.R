# TICKET-039 (beezdemand tracker): brms tier formula builders -------------------
#
# Verifies .dd_brms_formula() / .dd_brms_choice_formula() against the design
# spec (beezdemand internal_docs/design/DESIGN-brms-tier.md sections 2.6-2.7).
# Structural Stan-code assertions via brms::make_stancode()/make_standata()
# (pure R, no toolchain).

skip_if_not_installed("brms")

squish <- function(x) gsub("[[:space:]]+", "", paste(deparse(x, width.cutoff = 500), collapse = ""))
mu_formula <- function(spec) spec$formula$formula

dd_test_data <- function(n_id = 6) {
  set.seed(11)
  delays <- c(1, 7, 30, 90, 180, 365)
  d <- expand.grid(id = factor(seq_len(n_id)), x = delays)
  k_i <- exp(log(0.02) + rnorm(n_id, 0, 0.4))
  mu <- 1 / (1 + k_i[d$id] * d$x)
  d$y <- pmin(pmax(mu + rnorm(nrow(d), 0, 0.05), 0.01), 0.99)
  d
}

choice_test_data <- function(n_id = 6) {
  set.seed(12)
  d <- expand.grid(id = factor(seq_len(n_id)), delay = c(1, 7, 30, 90), rep = 1:3)
  d$rel <- 2
  d$choice <- rbinom(nrow(d), 1, 0.5)
  d
}

# ------------------------------------------------------------------------------
# Indifference-point equations (design 2.6)
# ------------------------------------------------------------------------------

test_that("mazur beta formula carries the sltb-analog squish and identity link", {
  spec <- beezdiscounting:::.dd_brms_formula(equation = "mazur", family = "beta")

  expect_s3_class(spec$formula, "brmsformula")
  expect_identical(
    squish(mu_formula(spec)),
    "y~1/(10^6)+(1-2/(10^6))*(1/(1+exp(logk)*x))"
  )
  expect_identical(spec$family$family, "beta")
  expect_identical(spec$family$link, "identity")
  expect_identical(spec$nlpars, "logk")
  expect_false(spec$has_s)
})

test_that("exponential / green-myerson / rachlin mu kernels match TMB", {
  sp_exp <- beezdiscounting:::.dd_brms_formula(equation = "exponential", family = "beta")
  expect_identical(
    squish(mu_formula(sp_exp)),
    "y~1/(10^6)+(1-2/(10^6))*(exp(-exp(logk)*x))"
  )

  sp_gm <- beezdiscounting:::.dd_brms_formula(equation = "green-myerson", family = "beta")
  expect_identical(
    squish(mu_formula(sp_gm)),
    "y~1/(10^6)+(1-2/(10^6))*((1+exp(logk)*x)^(-exp(logs)))"
  )
  expect_setequal(sp_gm$nlpars, c("logk", "logs"))
  expect_true(sp_gm$has_s)

  # rachlin guards pow(0, s): xzero/xsafe data columns (Codex design fold)
  sp_ra <- beezdiscounting:::.dd_brms_formula(equation = "rachlin", family = "beta")
  expect_identical(
    squish(mu_formula(sp_ra)),
    "y~1/(10^6)+(1-2/(10^6))*(xzero+(1-xzero)/(1+exp(logk)*xsafe^exp(logs)))"
  )
  expect_setequal(sp_ra$derived_cols, c("xzero", "xsafe"))
})

test_that("gaussian family is exact TMB parity: raw mu, no squish", {
  spec <- beezdiscounting:::.dd_brms_formula(equation = "mazur", family = "gaussian")
  expect_identical(squish(mu_formula(spec)), "y~1/(1+exp(logk)*x)")
  expect_identical(spec$family$family, "gaussian")
})

test_that("zoib boundary swaps in zero_one_inflated_beta", {
  spec <- beezdiscounting:::.dd_brms_formula(
    equation = "mazur", family = "beta", boundary = "zoib"
  )
  expect_identical(spec$family$family, "zero_one_inflated_beta")
  expect_identical(spec$family$link, "identity")
})

test_that("sltb errors with a pointer to the beta analog", {
  expect_error(
    beezdiscounting:::.dd_brms_formula(equation = "mazur", family = "sltb"),
    "beta"
  )
})

test_that("IP formulas generate valid Stan code with brms parameter names", {
  d <- dd_test_data()
  for (eq in c("mazur", "exponential", "green-myerson", "rachlin")) {
    for (fam in c("beta", "gaussian")) {
      spec <- beezdiscounting:::.dd_brms_formula(equation = eq, family = fam)
      dd <- d
      if (eq == "rachlin") {
        dd$xzero <- as.numeric(dd$x == 0)
        dd$xsafe <- ifelse(dd$x == 0, 1, dd$x)
      }
      scode <- brms::make_stancode(spec$formula, data = dd, family = spec$family)
      expect_gt(nchar(scode), 1000)
      expect_match(scode, "vector\\[K_logk\\] b_logk")
      expect_identical(grepl("b_logs", scode), spec$has_s, info = paste(eq, fam))
      if (fam == "beta") {
        expect_match(scode, "real<lower=0> phi")
      }
    }
  }
})

test_that("logk carries the subject random effect; logs is population-level", {
  d <- dd_test_data()
  spec <- beezdiscounting:::.dd_brms_formula(equation = "green-myerson", family = "beta")
  expect_match(squish(spec$formula$pforms$logk), "\\(1\\|id\\)")
  expect_false(grepl("\\|", squish(spec$formula$pforms$logs)))
  sdat <- brms::make_standata(spec$formula, data = d, family = spec$family)
  expect_identical(as.integer(sdat$M_1), 1L)
})

# ------------------------------------------------------------------------------
# Choice (structural) model (design 2.7)
# ------------------------------------------------------------------------------

test_that("choice structural formula is the TMB logit with bernoulli family", {
  spec <- beezdiscounting:::.dd_brms_choice_formula(equation = "mazur")
  expect_identical(
    squish(mu_formula(spec)),
    "choice~exp(loggamma)*(rel*(1/(1+exp(logk)*delay))-1)"
  )
  expect_identical(spec$family$family, "bernoulli")
  expect_setequal(spec$nlpars, c("logk", "loggamma"))

  spec_int <- beezdiscounting:::.dd_brms_choice_formula(
    equation = "exponential", intercept = TRUE
  )
  expect_identical(
    squish(mu_formula(spec_int)),
    "choice~b0+exp(loggamma)*(rel*(exp(-exp(logk)*delay))-1)"
  )
  expect_true("b0" %in% spec_int$nlpars)
})

test_that("choice formula generates valid Stan code", {
  d <- choice_test_data()
  spec <- beezdiscounting:::.dd_brms_choice_formula(equation = "mazur")
  scode <- brms::make_stancode(spec$formula, data = d)
  expect_gt(nchar(scode), 1000)
  expect_match(scode, "b_loggamma")
})

test_that("choice formula carries the factor design on logk (TICKET-048)", {
  d <- choice_test_data()
  d$group <- factor(ifelse(as.integer(d$id) <= 3, "ctrl", "treat"))

  spec <- beezdiscounting:::.dd_brms_choice_formula(
    equation = "mazur", factors = "group", data = d
  )
  expect_match(squish(spec$formula$pforms$logk), "group")
  expect_match(squish(spec$formula$pforms$logk), "\\(1\\|id\\)")
  # gamma (and b0) stay population-level: the design is on logk only
  expect_false(grepl("\\|", squish(spec$formula$pforms$loggamma)))

  sdat <- brms::make_standata(spec$formula, data = d)
  expect_identical(as.integer(sdat$K_logk), 2L)

  # continuous covariate route shares build_fixed_rhs()
  d$age <- as.numeric(d$id) * 3 + 20
  spec_cov <- beezdiscounting:::.dd_brms_choice_formula(
    equation = "mazur", continuous_covariates = "age", data = d
  )
  sdat_cov <- brms::make_standata(spec_cov$formula, data = d)
  expect_identical(as.integer(sdat_cov$K_logk), 2L)
})

# ------------------------------------------------------------------------------
# Snapshot: pin the canonical generated code against brms upgrades
# ------------------------------------------------------------------------------

test_that("canonical dd stancode snapshot is stable", {
  d <- dd_test_data()
  spec <- beezdiscounting:::.dd_brms_formula(equation = "mazur", family = "beta")
  expect_snapshot(cat(brms::make_stancode(spec$formula, data = d, family = spec$family)))
})
