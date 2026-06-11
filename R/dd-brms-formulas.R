# brms tier: formula builders -----------------------------------------------------
#
# Translates the TMB discounting models (src/MixedDiscounting.h,
# src/ChoiceDiscounting.h) into brms nonlinear formulas. Design spec:
# beezdemand internal_docs/design/DESIGN-brms-tier.md sections 2.6-2.7.
# Convention: infrastructure ported from beezdemand's brms tier with the
# .dd_ prefix.

#' Gate brms availability for the Bayesian tier
#' @noRd
.dd_brms_check_installed <- function() {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for the Bayesian (brms) modeling tier.\n",
      "Install it with install.packages(\"brms\").",
      call. = FALSE
    )
  }
}

#' Build the brms nonlinear formula for an indifference-point model
#'
#' The mean functions match `fit_dd_tmb()` exactly (`k = exp(logk)`,
#' `s = exp(logs)`), with `logk` carrying the subject random intercept
#' (`k ~ 1`, the v1 scope) and `logs` population-level. For
#' `family = "beta"` the mean is linearly squished into
#' `(1e-6, 1 - 1e-6)` -- the differentiable analog of the TMB sltb clamp --
#' and brms's `Beta(link = "identity")` is used (mu is naturally in (0,1)
#' for k > 0, so the identity link has no rejection region). For
#' `family = "gaussian"` the raw mean function matches
#' `fit_dd_tmb(family = "gaussian")` wherever the TMB template's mu clamp
#' into `[1e-6, 1 - 1e-6]` does not bind -- i.e. everywhere except extreme
#' decay underflow (`mu` below 1e-6, e.g. `k * delay > ~14` for the
#' exponential equation); the brms mean is left unclamped because a hard
#' clamp is not differentiably expressible in the shared R/Stan formula
#' language. The Rachlin equation
#' guards `pow(0, s)` (zero delays are allowed by `.dd_validate_ip()`)
#' via precomputed `xzero`/`xsafe` data columns, mirroring the TMB
#' CondExpGt guard.
#'
#' @param equation One of "mazur", "exponential", "green-myerson", "rachlin".
#' @param family "beta" (default) or "gaussian". "sltb" errors with a
#'   pointer: brms has no scale-location-truncated beta; beta + squeeze is
#'   the analog.
#' @param boundary "squeeze" (default; Smithson-Verkuilen, applied by the
#'   fitter) or "zoib" (swaps in `zero_one_inflated_beta`, which changes
#'   the estimand: k then describes interior responses only).
#' @param factors,factor_interaction,continuous_covariates Fixed-effect
#'   design on `logk` (as in `fit_dd_tmb()`); `logs` stays population-level.
#' @param data Optional data frame for single-level-factor dropping.
#' @return list(formula, family, nlpars, has_s, response_var, derived_cols,
#'   equation, boundary).
#' @noRd
.dd_brms_formula <- function(
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  family = c("beta", "gaussian"),
  boundary = c("squeeze", "zoib", "error"),
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  data = NULL
) {
  .dd_brms_check_installed()
  equation <- match.arg(equation)
  if (identical(family, "sltb")) {
    stop(
      "brms has no scale-location-truncated beta (sltb) family. ",
      "Use family = \"beta\" (with the default boundary = \"squeeze\"), ",
      "the closest analog, or family = \"gaussian\" for exact TMB parity.",
      call. = FALSE
    )
  }
  family <- match.arg(family)
  boundary <- match.arg(boundary)

  has_s <- equation %in% c("green-myerson", "rachlin")
  mu_str <- switch(equation,
    "mazur" = "1 / (1 + exp(logk) * x)",
    "exponential" = "exp(-exp(logk) * x)",
    "green-myerson" = "(1 + exp(logk) * x)^(-exp(logs))",
    "rachlin" = "xzero + (1 - xzero) / (1 + exp(logk) * xsafe^exp(logs))"
  )

  if (family == "beta") {
    # Linear squish into (1e-6, 1 - 1e-6): the differentiable analog of the
    # TMB sltb mu clamp; also makes x = 0 rows (mu = 1) non-fatal.
    # Written as 1/(10^6) NOT 1e-6: brms's Stan-code printer inserts spaces
    # around '-' inside scientific literals ("1e - 06"), a stanc syntax
    # error that make_stancode() does not catch (it never parses).
    mu_str <- paste0("1/(10^6) + (1 - 2/(10^6)) * (", mu_str, ")")
    fam <- if (boundary == "zoib") {
      brms::zero_one_inflated_beta(link = "identity")
    } else {
      brms::Beta(link = "identity")
    }
  } else {
    fam <- brms::brmsfamily("gaussian", link = "identity")
  }

  # build_fixed_rhs() returns a formula here; deparse before string assembly
  # (sub() on a formula coerces to a length-2 character and triggers the
  # deprecated formula.character path downstream).
  rhs <- build_fixed_rhs(
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    data = data
  )
  rhs_core <- sub(
    "^~\\s*", "",
    paste(deparse(stats::as.formula(rhs), width.cutoff = 500), collapse = " ")
  )
  pforms <- list(stats::as.formula(paste0("logk ~ ", rhs_core, " + (1 | id)")))
  if (has_s) {
    pforms <- c(pforms, list(stats::as.formula("logs ~ 1")))
  }

  bform <- do.call(
    brms::bf,
    c(list(stats::as.formula(paste("y ~", mu_str))), pforms, list(nl = TRUE))
  )

  list(
    formula = bform,
    family = fam,
    nlpars = c("logk", if (has_s) "logs"),
    has_s = has_s,
    response_var = "y",
    derived_cols = if (equation == "rachlin") c("xzero", "xsafe") else character(0),
    equation = equation,
    boundary = if (family == "beta") boundary else NA_character_
  )
}

#' Build the brms formula for the structural choice model
#'
#' The TMB structural likelihood (src/ChoiceDiscounting.h):
#' `logit P(LL) = [b0] + gamma * ((ll/ss) * D(k, delay) - 1)` with
#' `gamma = exp(loggamma)` and `k = exp(logk)` carrying the subject random
#' intercept. With `bernoulli("logit")` the nonlinear formula IS the logit,
#' matching TMB exactly. The fitter precomputes `rel = ll/ss`.
#'
#' @param equation "mazur" or "exponential".
#' @param intercept Include the optional bias term `b0`.
#' @return list(formula, family, nlpars, response_var, derived_cols,
#'   equation, intercept).
#' @noRd
.dd_brms_choice_formula <- function(
  equation = c("mazur", "exponential"),
  intercept = FALSE
) {
  .dd_brms_check_installed()
  equation <- match.arg(equation)

  disc_str <- switch(equation,
    "mazur" = "1 / (1 + exp(logk) * delay)",
    "exponential" = "exp(-exp(logk) * delay)"
  )
  rhs <- paste0("exp(loggamma) * (rel * (", disc_str, ") - 1)")
  if (isTRUE(intercept)) {
    rhs <- paste0("b0 + ", rhs)
  }

  pforms <- list(
    stats::as.formula("logk ~ 1 + (1 | id)"),
    stats::as.formula("loggamma ~ 1")
  )
  if (isTRUE(intercept)) {
    pforms <- c(pforms, list(stats::as.formula("b0 ~ 1")))
  }

  bform <- do.call(
    brms::bf,
    c(list(stats::as.formula(paste("choice ~", rhs))), pforms, list(nl = TRUE))
  )

  list(
    formula = bform,
    family = brms::bernoulli(link = "logit"),
    nlpars = c("logk", "loggamma", if (isTRUE(intercept)) "b0"),
    response_var = "choice",
    derived_cols = "rel",
    equation = equation,
    intercept = isTRUE(intercept)
  )
}
