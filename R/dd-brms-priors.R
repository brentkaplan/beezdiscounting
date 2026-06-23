# brms tier: default priors --------------------------------------------------------
#
# Weakly-informative defaults on the estimation (natural-log) scale. Design:
# beezdemand internal_docs/design/DESIGN-brms-tier.md section 3.2.

#' Contract number format for autoscaled prior strings
#' @noRd
.dd_brms_fmt_num <- function(x) {
  format(x, digits = 6, scientific = FALSE, trim = TRUE)
}

#' Default priors for Bayesian (brms) delay-discounting models
#'
#' Returns the default prior table used by `fit_dd_brms()`, so the defaults
#' can be inspected, modified row-wise, and passed back via the fitter's
#' `prior` argument. The `logk` location is the principled fix for k's
#' delay-unit dependence: with `autoscale = TRUE` (the default whenever
#' `data` is supplied) it centers `k * median(delay) = 1` -- the delay at
#' which the Mazur curve crosses 0.5 -- via `normal(-log(median(x)), 2.5)`;
#' the static fallback is `normal(-4.5, 2.5)`. The anchors used are
#' attached as `attr(, "autoscale_info")`; numeric values are formatted
#' with `format(x, digits = 6, scientific = FALSE)`.
#'
#' Other defaults: `logs ~ normal(0, 0.5)` (two-parameter equations; s is
#' near 1 a priori), `sd(logk) ~ student_t(3, 0, 1)`,
#' `phi ~ gamma(2, 0.1)` (Beta precision; mean 20, far from brms's
#' near-improper `gamma(0.01, 0.01)`), and
#' `sigma ~ student_t(3, 0, 0.25)` for the Gaussian family (y is a
#' proportion in the unit interval, so sd(y) <= 0.5).
#'
#' @param equation Discounting equation (TMB-tier vocabulary).
#' @param family `"beta"` or `"gaussian"`.
#' @param data Optional data frame used for autoscaling.
#' @param y_var,x_var Column names in `data` (canonical defaults).
#' @param factors,factor_interaction,continuous_covariates Fixed-effect
#'   design on `logk`, as passed to the fitter. When the design has
#'   non-intercept coefficients (derived through `build_fixed_rhs()`, so
#'   single-level dropped factors do not count), a fold-change
#'   `normal(0, 1)` class-level coefficient prior is added; with an
#'   intercept-only design it is omitted (it would be unused, and brms
#'   warns).
#' @param autoscale Logical; defaults to `TRUE` when `data` is supplied.
#' @param random_effects `k ~ 1` (default) or `k + phi ~ 1`. The latter
#'   re-keys the beta precision: phi becomes a predicted distributional
#'   parameter, so the scalar `gamma(2, 0.1)` is replaced by a log-scale
#'   intercept prior, a half-t precision-RE SD, and (for `"pdSymm"`) an LKJ
#'   correlation prior.
#' @param covariance_structure `(log k, log phi)` covariance for `k + phi ~ 1`:
#'   `"pdSymm"` (default, correlated) or `"pdDiag"` (independent).
#'
#' @return A `brmsprior` data frame, with `attr(, "autoscale_info")` when
#'   autoscaling was used.
#' @export
default_dd_priors <- function(
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  family = c("beta", "gaussian"),
  data = NULL,
  y_var = "y",
  x_var = "x",
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  autoscale = !is.null(data),
  random_effects = k ~ 1,
  covariance_structure = c("pdSymm", "pdDiag")
) {
  .dd_brms_check_installed()
  equation <- match.arg(equation)
  family <- match.arg(family)
  re <- .dd_normalize_re(
    random_effects,
    covariance_structure = covariance_structure
  )
  phi_re <- length(re$blocks[[1]]$param) == 2L
  correlated <- phi_re && identical(re$blocks[[1]]$pdmat_class, "pdSymm")
  if (phi_re && family != "beta") {
    stop(
      "phi random effects (`k + phi ~ 1`) require family = \"beta\" (the ",
      "gaussian family has a residual SD, not a precision).",
      call. = FALSE
    )
  }
  fmt <- .dd_brms_fmt_num
  has_s <- equation %in% c("green-myerson", "rachlin")

  # build_fixed_rhs() returns a FORMULA in this package (beezdemand's port
  # returns a string): count non-intercept terms rather than comparing
  # representations.
  rhs <- build_fixed_rhs(
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    data = data
  )
  has_coefs <- length(attr(
    stats::terms(stats::as.formula(rhs)),
    "term.labels"
  )) >
    0

  info <- NULL
  if (isTRUE(autoscale)) {
    if (is.null(data)) {
      stop("`autoscale` = TRUE requires `data`.", call. = FALSE)
    }
    x <- data[[x_var]]
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0L) {
      stop("Cannot autoscale priors: no positive delays.", call. = FALSE)
    }
    info <- list(median_delay = stats::median(x))
  }

  p_logk <- if (!is.null(info)) {
    paste0("normal(", fmt(-log(info$median_delay)), ", 2.5)")
  } else {
    "normal(-4.5, 2.5)"
  }

  rows <- list(
    brms::set_prior(p_logk, class = "b", coef = "Intercept", nlpar = "logk"),
    brms::set_prior("student_t(3, 0, 1)", class = "sd", nlpar = "logk")
  )
  if (has_coefs) {
    rows <- c(
      rows,
      list(
        brms::set_prior("normal(0, 1)", class = "b", nlpar = "logk")
      )
    )
  }
  if (has_s) {
    rows <- c(
      rows,
      list(
        brms::set_prior(
          "normal(0, 0.5)",
          class = "b",
          coef = "Intercept",
          nlpar = "logs"
        )
      )
    )
  }
  if (family == "beta") {
    if (phi_re) {
      # phi is predicted on the log link: a location prior on the (log)
      # intercept (centered at log of the scalar default's prior mean, 20), a
      # half-t SD on the precision RE, and -- when correlated -- an LKJ prior on
      # the (log k, log phi) correlation.
      rows <- c(
        rows,
        list(
          brms::set_prior(
            paste0("student_t(3, ", fmt(log(20)), ", 1)"),
            class = "Intercept",
            dpar = "phi"
          ),
          brms::set_prior("student_t(3, 0, 1)", class = "sd", dpar = "phi")
        )
      )
      if (correlated) {
        rows <- c(
          rows,
          list(
            brms::set_prior("lkj(2)", class = "cor", group = "id")
          )
        )
      }
    } else {
      rows <- c(rows, list(brms::set_prior("gamma(2, 0.1)", class = "phi")))
    }
  } else {
    rows <- c(
      rows,
      list(brms::set_prior("student_t(3, 0, 0.25)", class = "sigma"))
    )
  }

  pri <- do.call(c, rows)
  if (!is.null(info)) {
    attr(pri, "autoscale_info") <- info
  }
  pri
}

#' Default priors for the Bayesian (brms) choice model
#'
#' `loggamma ~ normal(1, 1)` (choice sensitivity, typically 1-20) and the
#' optional logit-scale bias `b0 ~ normal(0, 1.5)`; `logk` defaults as in
#' [default_dd_priors()] (autoscaled to the median delay when `data` is
#' supplied; `delay_var` names the delay column).
#'
#' @param equation `"mazur"` or `"exponential"`.
#' @param intercept Include the `b0` prior.
#' @param data Optional data frame used for `logk` autoscaling.
#' @param delay_var Delay column name in `data`.
#' @param factors,factor_interaction,continuous_covariates Fixed-effect
#'   design on `logk`, as passed to [fit_dd_choice_brms()]. When the design
#'   has non-intercept coefficients (derived through `build_fixed_rhs()`),
#'   a fold-change `normal(0, 1)` class-level coefficient prior is added;
#'   with an intercept-only design it is omitted (it would be unused, and
#'   brms warns).
#' @param autoscale Logical; defaults to `TRUE` when `data` is supplied.
#' @return A `brmsprior` data frame.
#' @export
default_dd_choice_priors <- function(
  equation = c("mazur", "exponential"),
  intercept = FALSE,
  data = NULL,
  delay_var = "delay",
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  autoscale = !is.null(data)
) {
  .dd_brms_check_installed()
  equation <- match.arg(equation)
  fmt <- .dd_brms_fmt_num

  # Count non-intercept terms through the same build_fixed_rhs() route as the
  # formula builder (single-level dropped factors do not count).
  fe_rhs <- build_fixed_rhs(
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    data = data
  )
  has_coefs <-
    length(attr(stats::terms(stats::as.formula(fe_rhs)), "term.labels")) > 0

  info <- NULL
  if (isTRUE(autoscale)) {
    if (is.null(data)) {
      stop("`autoscale` = TRUE requires `data`.", call. = FALSE)
    }
    x <- data[[delay_var]]
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0L) {
      stop("Cannot autoscale priors: no positive delays.", call. = FALSE)
    }
    info <- list(median_delay = stats::median(x))
  }

  p_logk <- if (!is.null(info)) {
    paste0("normal(", fmt(-log(info$median_delay)), ", 2.5)")
  } else {
    "normal(-4.5, 2.5)"
  }

  rows <- list(
    brms::set_prior(p_logk, class = "b", coef = "Intercept", nlpar = "logk"),
    brms::set_prior("student_t(3, 0, 1)", class = "sd", nlpar = "logk"),
    brms::set_prior(
      "normal(1, 1)",
      class = "b",
      coef = "Intercept",
      nlpar = "loggamma"
    )
  )
  if (has_coefs) {
    rows <- c(
      rows,
      list(
        brms::set_prior("normal(0, 1)", class = "b", nlpar = "logk")
      )
    )
  }
  if (isTRUE(intercept)) {
    rows <- c(
      rows,
      list(
        brms::set_prior(
          "normal(0, 1.5)",
          class = "b",
          coef = "Intercept",
          nlpar = "b0"
        )
      )
    )
  }

  pri <- do.call(c, rows)
  if (!is.null(info)) {
    attr(pri, "autoscale_info") <- info
  }
  pri
}

#' Merge user priors over defaults (user rows win on the full prior key)
#'
#' Verbatim port of beezdemand's `.brms_merge_priors()` (the .dd_ convention).
#' @noRd
.dd_brms_merge_priors <- function(user, defaults) {
  if (is.null(user)) {
    return(defaults)
  }
  if (!inherits(user, "brmsprior")) {
    stop(
      "`prior` must be a brmsprior object (see brms::set_prior()).",
      call. = FALSE
    )
  }
  key <- function(p) {
    paste(p$class, p$coef, p$group, p$dpar, p$nlpar, p$resp, sep = "\r")
  }
  keep <- !(key(defaults) %in% key(user))
  out <- rbind(user, defaults[keep, , drop = FALSE])
  rownames(out) <- NULL
  attr(out, "autoscale_info") <- attr(defaults, "autoscale_info")
  out
}
