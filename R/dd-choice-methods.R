# ==============================================================================
# S3 Methods for beezdiscounting_choice Objects (structural SS-vs-LL GLMM)
#
# Mirrors R/dd-tmb-methods.R (the IP family's S3 surface). The choice coefficient
# vector orders `beta_k` (one per design column) FIRST, then `log_sigma_u`,
# `log_gamma`, and `beta0` (only when intercept = TRUE). Display term mapping:
#   beta_k       -> k:<design colname>   (so .dd_transform_coef_table back-
#                                         transforms it via its ^k($|_|:) regex)
#   log_gamma    -> gamma                (choice sensitivity; transformed
#                                         EXPLICITLY via .dd_transform_est_se)
#   beta0        -> beta0                (identity scale; NEVER transformed)
#   log_sigma_u  -> log_sigma_u          (variance component)
# ==============================================================================

#' Build display term names from a beezdiscounting_choice coefficient vector
#'
#' Maps raw optimizer names to readable display names. `beta_k` columns become
#' `k:<design colname>` using `formula_details$X` colnames; `log_gamma` becomes
#' `gamma`; every other coefficient (`log_sigma_u`, `beta0`) keeps its raw name.
#' Shared by `tidy()`, `summary()`, `confint()`, `print()`.
#'
#' @param object A `beezdiscounting_choice` object.
#' @param nms Character vector of raw parameter names (default
#'   `names(object$model$coefficients)`).
#' @return List with `term` (display names), `k_idx` (beta_k positions),
#'   `other_idx` (non-beta positions).
#' @keywords internal
.dd_choice_term_names <- function(object, nms = NULL) {
  if (is.null(nms)) {
    nms <- names(object$model$coefficients)
  }
  k_idx <- which(nms == "beta_k")
  other_idx <- which(nms != "beta_k")

  k_colnames <- colnames(object$formula_details$X)
  if (is.null(k_colnames)) {
    k_colnames <- paste0("X", seq_along(k_idx))
  }

  term <- character(length(nms))
  term[k_idx] <- paste0("k:", k_colnames)
  term[other_idx] <- nms[other_idx]
  # Display log_gamma as "gamma" for the coefficient/summary tables. The
  # param-space transformer does NOT key on gamma (only k/s/phi), so gamma is
  # back-transformed EXPLICITLY by the callers (R6).
  term[nms == "log_gamma"] <- "gamma"

  list(
    term = term,
    k_idx = k_idx,
    other_idx = as.integer(other_idx)
  )
}


#' Named standard-error vector aligned to the coefficient vector
#'
#' Returns `object$model$se` aligned to `names(object$model$coefficients)`.
#' Returns `NA` SEs when `se_available` is `FALSE` (non-PD Hessian or failed
#' sdreport), so confint()/tidy()/summary() never present unreliable Wald
#' intervals/p-values.
#'
#' @param object A `beezdiscounting_choice` object.
#' @return Named numeric vector parallel to `object$model$coefficients`.
#' @keywords internal
.dd_choice_model_se <- function(object) {
  co <- object$model$coefficients
  na_se <- stats::setNames(rep(NA_real_, length(co)), names(co))
  if (isFALSE(object$se_available)) {
    return(na_se)
  }
  se <- object$model$se
  if (is.null(se)) {
    return(na_se)
  }
  # Defensive: align to coefficient order/names. The choice se vector is built
  # parallel to opt$par (same names) in .dd_choice_extract_estimates, so this is
  # a positional/name identity in practice.
  stats::setNames(unname(se)[match(names(co), names(se))], names(co))
}


# --- logLik / AIC / BIC / nobs ---

#' @importFrom stats AIC BIC logLik nobs
#' @export
logLik.beezdiscounting_choice <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- length(object$opt$par)
  attr(ll, "nobs") <- object$param_info$n_obs
  class(ll) <- "logLik"
  ll
}

#' @export
AIC.beezdiscounting_choice <- function(object, ..., k = 2) {
  if (k != 2) {
    nll <- -object$loglik
    n_params <- length(object$opt$par)
    return(2 * nll + k * n_params)
  }
  object$AIC
}

#' @export
BIC.beezdiscounting_choice <- function(object, ...) {
  object$BIC
}

#' @export
nobs.beezdiscounting_choice <- function(object, ...) {
  object$param_info$n_obs %||% nrow(object$data)
}


# --- coef / fixef / ranef ---

#' Extract coefficients from a structural choice model
#'
#' Returns the optimizer's flat named parameter vector: `beta_k` (one entry per
#' fixed-effect design column, on the log-k scale), `log_sigma_u`, `log_gamma`,
#' and `beta0` (only when `intercept = TRUE`). This is the numeric escape hatch
#' consumed by tooling such as `car::deltaMethod`.
#'
#' @param object A `beezdiscounting_choice` object.
#' @param ... Unused.
#' @return Named numeric vector.
#' @importFrom stats coef
#' @export
coef.beezdiscounting_choice <- function(object, ...) {
  object$model$coefficients
}

#' Extract fixed effects from a structural choice model
#'
#' @param object A `beezdiscounting_choice` object.
#' @param ... Unused.
#' @return Named numeric vector (identical to `coef()`).
#' @importFrom nlme fixef
#' @export
fixef.beezdiscounting_choice <- function(object, ...) {
  coef(object)
}

#' Extract subject-level random effects from a structural choice model
#'
#' @param object A `beezdiscounting_choice` object.
#' @param ... Unused.
#' @return Data frame with `id`, the standardized random-intercept deviate
#'   `u_i` (such that `log k_i = X beta + sigma_u * u_i`), and the resolved
#'   per-subject discount rate `k`.
#' @importFrom nlme ranef
#' @export
ranef.beezdiscounting_choice <- function(object, ...) {
  sp <- object$subject_pars
  keep <- intersect(c("id", "u_i", "k"), names(sp))
  out <- sp[, keep, drop = FALSE]
  rownames(out) <- NULL
  out
}


# --- predict ---

#' Predict from a structural choice discounting model
#'
#' @param object A `beezdiscounting_choice` object.
#' @param newdata Optional data frame. `NULL` uses the fitting data. A supplied
#'   `newdata` must use the package's canonical column names (`id`, `ss_amount`,
#'   `ll_amount`, `delay`, plus any factor/covariate columns).
#' @param type `"prob"` (default; fitted choice probabilities P(LL) on `(0,1)`)
#'   or `"parameters"` (the per-subject parameter tibble).
#' @param level For `type = "prob"`: `"subject"` (default; conditions on each
#'   subject's estimated random intercept, requires the id column) or
#'   `"population"` (random effects set to zero - the population-mean curve; no
#'   id column needed).
#' @param ... Unused.
#'
#' @return
#' - `type = "parameters"`: the subject-parameter tibble (one row per subject).
#' - `type = "prob"`: `newdata` as a tibble plus a `.prob` column of fitted
#'   choice probabilities.
#'
#' @export
predict.beezdiscounting_choice <- function(object,
                                           newdata = NULL,
                                           type = c("prob", "parameters"),
                                           level = "subject",
                                           ...) {
  type <- match.arg(type)
  if (!is.character(level)) {
    cli::cli_abort(
      c("{.arg level} must be a character vector.",
        i = "{.arg level} should be one of {.val \"subject\"} or {.val \"population\"}."),
      call = NULL
    )
  }
  level <- match.arg(level, c("subject", "population"))

  if (type == "parameters") {
    return(tibble::as_tibble(object$subject_pars))
  }

  # type == "prob"
  equation <- object$param_info$equation
  coefs    <- object$model$coefficients
  gamma    <- exp(unname(coefs[["log_gamma"]]))
  beta0    <- if (isTRUE(object$param_info$intercept)) unname(coefs[["beta0"]]) else 0
  if (is.null(newdata)) newdata <- object$data
  out <- tibble::as_tibble(newdata)

  # Supplied newdata uses the CANONICAL trial columns. Fail cleanly (not with an
  # opaque length-0 result) if any structural column is absent. Mirrors the IP
  # predict guard.
  needed_trial <- c("ss_amount", "ll_amount", "delay")
  missing_trial <- setdiff(needed_trial, names(newdata))
  if (length(missing_trial) > 0L) {
    cli::cli_abort(c(
      "{.arg newdata} is missing trial column{?s}: {.val {missing_trial}}.",
      "i" = "Supplied {.arg newdata} uses the canonical trial names \\
             {.val ss_amount} / {.val ll_amount} / {.val delay}."
    ))
  }

  k_row <- .dd_tmb_predict_k(object, newdata, level = level)
  eta <- .dd_choice_structural_eta(
    k = k_row, ss_amount = newdata$ss_amount, ll_amount = newdata$ll_amount,
    delay = newdata$delay, equation = equation, gamma = gamma, beta0 = beta0
  )
  out$.prob <- stats::plogis(eta)
  out
}


# --- fitted / residuals / augment ---

#' Shared fitted-probabilities / residuals back end
#'
#' Calls `predict(type = "prob")` at the requested level and returns the fitted
#' choice probabilities, raw residuals (`choice - p`), and the data frame used.
#'
#' @param object A `beezdiscounting_choice` fit.
#' @param newdata Optional data frame; `NULL` -> training data.
#' @param level `"subject"` (default) or `"population"`.
#' @return List with `.fitted`, `.resid`, and `data`.
#' @keywords internal
.dd_choice_fitted_resid <- function(object, newdata = NULL,
                                    level = c("subject", "population")) {
  level     <- match.arg(level)
  data_used <- if (is.null(newdata)) object$data else newdata
  pred      <- predict(object, newdata = data_used, type = "prob", level = level)
  fitted_vals <- pred$.prob
  y_obs <- data_used[[object$param_info$y_var]]
  if (is.null(y_obs)) {
    cli::cli_abort(c(
      "{.arg newdata} must contain the response column \\
       {.val {object$param_info$y_var}} to compute residuals.",
      "i" = "Supplied {.arg newdata} uses the canonical name {.val choice}."
    ))
  }
  list(.fitted = fitted_vals, .resid = y_obs - fitted_vals, data = data_used)
}


#' Fitted values for a beezdiscounting_choice fit
#'
#' @param object A `beezdiscounting_choice` object.
#' @param level `"subject"` (default) or `"population"`.
#' @param ... Unused.
#' @return Numeric vector of fitted choice probabilities, length `nobs(object)`.
#' @export
fitted.beezdiscounting_choice <- function(object,
                                          level = c("subject", "population"),
                                          ...) {
  level <- match.arg(level)
  .dd_choice_fitted_resid(object, level = level)$.fitted
}


#' Residuals for a beezdiscounting_choice fit
#'
#' @param object A `beezdiscounting_choice` object.
#' @param type `"response"` (default; raw `choice - p`) or `"pearson"`
#'   (divided by `sqrt(p * (1 - p))`).
#' @param level `"subject"` (default) or `"population"`.
#' @param ... Unused.
#' @return Numeric vector of residuals, length `nobs(object)`.
#' @export
residuals.beezdiscounting_choice <- function(object,
                                             type  = c("response", "pearson"),
                                             level = c("subject", "population"),
                                             ...) {
  type  <- match.arg(type)
  level <- match.arg(level)
  fr    <- .dd_choice_fitted_resid(object, level = level)
  if (type == "response") return(fr$.resid)
  # Clamp the probability in the Pearson SD denominator so a fitted prob that
  # saturates to 0/1 does not divide-by-zero. The response residual stays exact.
  p_sd <- pmin(pmax(fr$.fitted, 1e-6), 1 - 1e-6)
  fr$.resid / sqrt(p_sd * (1 - p_sd))
}


#' Augment a beezdiscounting_choice model
#'
#' Returns the data used for fitting (or `newdata`) as a tibble with three
#' diagnostic columns appended:
#' - `.fitted`: subject-conditional fitted choice probability P(LL).
#' - `.resid`: raw (response) residual `choice - .fitted`.
#' - `.std_resid`: Pearson residual `.resid / sqrt(p * (1 - p))`.
#'
#' @param x A `beezdiscounting_choice` object.
#' @param newdata Optional data frame. `NULL` uses the fitting data.
#' @param ... Unused.
#' @return A tibble with the same rows as the data plus `.fitted`, `.resid`,
#'   and `.std_resid`.
#'
#' @importFrom generics augment
#' @export
augment.beezdiscounting_choice <- function(x, newdata = NULL, ...) {
  fr  <- .dd_choice_fitted_resid(x, newdata = newdata, level = "subject")
  out <- tibble::as_tibble(fr$data)
  out$.fitted    <- fr$.fitted
  out$.resid     <- fr$.resid
  # Clamp the probability in the Pearson SD denominator so a fitted prob that
  # saturates to 0/1 does not divide-by-zero. The response residual stays exact.
  p_sd <- pmin(pmax(fr$.fitted, 1e-6), 1 - 1e-6)
  out$.std_resid <- fr$.resid / sqrt(p_sd * (1 - p_sd))
  out
}


# --- tidy / glance ---

#' Variance components for a beezdiscounting_choice fit
#'
#' Reports the random-intercept SD on the log10-k scale (divides the natural-log
#' SD by `log(10)` for comparability with `nlme::VarCorr()`).
#'
#' @param object A `beezdiscounting_choice` fit.
#' @return A data frame with columns `Component`, `Estimate`, `Scale`.
#' @keywords internal
.dd_choice_variance_components <- function(object) {
  coefs <- object$model$coefficients
  ln10  <- log(10)
  data.frame(
    Component = "sigma_u (log10-k RE SD)",
    Estimate  = exp(coefs[["log_sigma_u"]]) / ln10,
    Scale     = "log10",
    stringsAsFactors = FALSE
  )
}


#' Tidy a beezdiscounting_choice model into a coefficient tibble
#'
#' Returns fixed-effect rows (the log-k coefficients), the shape rows (the
#' choice-sensitivity `gamma` and, when present, the choice-bias `beta0`), and
#' the variance-component row, following the broom coefficient-table contract.
#'
#' `estimate` and `std.error` are reported on the `report_space` scale for the
#' fixed-effect (`beta_k`) rows and for `gamma` (which is transformed
#' EXPLICITLY since the param-space transformer keys only on k/s/phi). `beta0`
#' is on the identity (logit-intercept) scale and is NEVER transformed across
#' report spaces. `statistic` and `p.value` are always computed on the
#' estimation (internal) scale - Wald statistics are not recomputed after
#' back-transforming (broom convention).
#'
#' @param x A `beezdiscounting_choice` object.
#' @param effects Character vector: `"fixed"` (log-k fixed-effect rows + the
#'   shape rows gamma/beta0), `"ran_pars"` (the RE SD), or both (default).
#' @param report_space `"natural"`, `"log10"`, `"internal"`, or `"log"` -
#'   reporting scale for the fixed-effect `estimate`/`std.error`. Default is
#'   `"natural"`.
#' @param ... Unused.
#' @return A tibble with exactly 8 columns in this order: `term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, `component`, `estimate_scale`,
#'   `term_display`. Fixed-effect rows carry `component == "fixed"`; gamma/beta0
#'   carry `component == "shape"`; variance rows carry `component == "variance"`.
#'
#' @importFrom generics tidy
#' @export
tidy.beezdiscounting_choice <- function(x,
                                        effects      = c("fixed", "ran_pars"),
                                        report_space = c("natural", "log10", "internal", "log"),
                                        ...) {
  effects      <- match.arg(effects, several.ok = TRUE)
  report_space <- match.arg(report_space)

  result <- tibble::tibble(
    term           = character(),
    estimate       = numeric(),
    std.error      = numeric(),
    statistic      = numeric(),
    p.value        = numeric(),
    component      = character(),
    estimate_scale = character(),
    term_display   = character()
  )

  if ("fixed" %in% effects) {
    coefs <- x$model$coefficients
    se    <- .dd_choice_model_se(x)
    nms   <- names(coefs)
    tn    <- .dd_choice_term_names(x, nms)

    is_fixed <- nms == "beta_k"

    # Wald statistics on the estimation (internal) scale - never recomputed
    # after back-transforming (broom convention).
    z_val <- coefs / se
    p_val <- 2 * stats::pnorm(-abs(z_val))

    fixed <- tibble::tibble(
      term           = tn$term[is_fixed],
      estimate       = unname(coefs[is_fixed]),
      std.error      = unname(se[is_fixed]),
      statistic      = unname(z_val[is_fixed]),
      p.value        = unname(p_val[is_fixed]),
      component      = "fixed",
      estimate_scale = "log",     # internal space for beta_k is log-k
      term_display   = tn$term[is_fixed]
    )
    # Back-transform estimate + std.error to report_space (k rows only; the
    # transformer keys on ^k($|_|:)). statistic/p.value are left on the
    # estimation scale.
    fixed <- .dd_transform_coef_table(
      coef_tbl       = fixed,
      report_space   = report_space,
      internal_space = "log"
    )
    fixed <- fixed[, setdiff(names(fixed), "estimate_internal"), drop = FALSE]
    result <- dplyr::bind_rows(result, fixed)

    # gamma (log_gamma): choice-sensitivity shape parameter. R6 - transform
    # EXPLICITLY (the param-space transformer does not key on gamma).
    g_pos <- which(nms == "log_gamma")
    z_g   <- coefs[g_pos] / se[g_pos]
    p_g   <- 2 * stats::pnorm(-abs(z_g))
    to_g  <- if (report_space == "internal") "log" else report_space
    g_tr  <- .dd_transform_est_se(
      estimate = unname(coefs[g_pos]), se = unname(se[g_pos]),
      from = "log", to = to_g
    )
    # to_g is exhaustive over {natural, log10, log}; the stop() is a defensive
    # guard that can never fire given the upstream match.arg on report_space.
    g_disp <- switch(to_g,
      natural = "gamma",
      log10   = "log10(gamma)",
      log     = "log(gamma)",
      stop("unexpected report space '", to_g, "'", call. = FALSE)
    )
    shape <- tibble::tibble(
      term           = "gamma",
      estimate       = g_tr$estimate,
      std.error      = g_tr$se,
      statistic      = unname(z_g),
      p.value        = unname(p_g),
      component      = "shape",
      estimate_scale = to_g,
      term_display   = g_disp
    )
    result <- dplyr::bind_rows(result, shape)

    # beta0 (choice bias, only when intercept = TRUE). R6 - identity scale,
    # NEVER transformed/exponentiated across report spaces.
    if (isTRUE(x$param_info$intercept)) {
      b_pos <- which(nms == "beta0")
      z_b   <- coefs[b_pos] / se[b_pos]
      p_b   <- 2 * stats::pnorm(-abs(z_b))
      beta0_row <- tibble::tibble(
        term           = "beta0",
        estimate       = unname(coefs[b_pos]),
        std.error      = unname(se[b_pos]),
        statistic      = unname(z_b),
        p.value        = unname(p_b),
        component      = "shape",
        estimate_scale = "identity",
        term_display   = "beta0"
      )
      result <- dplyr::bind_rows(result, beta0_row)
    }
  }

  if ("ran_pars" %in% effects) {
    vc  <- .dd_choice_variance_components(x)
    ran <- tibble::tibble(
      term           = vc$Component,
      estimate       = vc$Estimate,
      std.error      = NA_real_,
      statistic      = NA_real_,
      p.value        = NA_real_,
      component      = "variance",
      estimate_scale = vc$Scale,
      term_display   = vc$Component
    )
    result <- dplyr::bind_rows(result, ran)
  }

  if (isFALSE(x$converged) || isFALSE(x$se_available)) {
    attr(result, "se_warning") <-
      "Standard errors / convergence are unreliable; CIs and p-values may be invalid."
  }

  result
}


#' Glance at a beezdiscounting_choice model
#'
#' Returns a one-row tibble of model-level summary statistics suitable for use
#' with `broom::glance()` workflows and multi-model comparison tables.
#'
#' @param x A `beezdiscounting_choice` object.
#' @param ... Unused.
#' @return A one-row tibble with columns `model_class`, `backend`, `mode`,
#'   `equation`, `nobs`, `n_subjects`, `n_random_effects`, `converged`,
#'   `logLik`, `AIC`, `BIC`.
#'
#' @importFrom generics glance
#' @export
glance.beezdiscounting_choice <- function(x, ...) {
  tibble::tibble(
    model_class      = "beezdiscounting_choice",
    backend          = "TMB_choice",
    mode             = x$param_info$mode,
    equation         = x$param_info$equation,
    nobs             = x$param_info$n_obs,
    n_subjects       = x$param_info$n_subjects,
    n_random_effects = x$param_info$n_random_effects,
    converged        = x$converged,
    logLik           = x$loglik,
    AIC              = x$AIC,
    BIC              = x$BIC
  )
}


# --- confint ---

#' Confidence intervals for a structural choice discounting model
#'
#' Wald (Hessian-based) confidence intervals: `estimate +/- z * se` on the
#' internal scale, then back-transformed to `report_space`.
#'
#' @param object A `beezdiscounting_choice` object.
#' @param parm Optional character vector for filtering. Accepts display names
#'   (`"k:(Intercept)"`, `"gamma"`) **or** raw optimizer names (`"beta_k"`,
#'   `"log_gamma"`, `"log_sigma_u"`, `"beta0"`). `NULL` returns all
#'   coefficients.
#' @param level Confidence level (default `0.95`).
#' @param report_space `"internal"` (default; all coefficients on their
#'   estimation scale) or `"natural"` (exponentiate the `beta_k` rows so the
#'   intercept is `k` at the reference level, and exponentiate `log_gamma` so
#'   the row is `gamma`; `log_sigma_u` and `beta0` stay on their internal
#'   scales).
#' @param ... Unused.
#' @return A tibble with columns `term`, `estimate`, `conf.low`, `conf.high`,
#'   `level`.
#'
#' @exportS3Method stats::confint beezdiscounting_choice
confint.beezdiscounting_choice <- function(object,
                                           parm         = NULL,
                                           level        = 0.95,
                                           report_space = c("internal", "natural"),
                                           ...) {
  report_space <- match.arg(report_space)
  coefs  <- object$model$coefficients
  se_vec <- .dd_choice_model_se(object)
  if (isFALSE(object$se_available)) {
    cli::cli_warn(c(
      "!" = "Standard errors are unreliable (non-PD Hessian or sdreport \\
             failed); confidence intervals are returned as {.val NA}.",
      "i" = "Check convergence and model identifiability before interpreting \\
             uncertainty."
    ))
  }
  nms  <- names(coefs)
  tn   <- .dd_choice_term_names(object, nms)
  term <- tn$term

  # Filter by display name OR raw name (dual-name matching).
  if (!is.null(parm)) {
    keep   <- term %in% parm | nms %in% parm
    coefs  <- coefs[keep]
    se_vec <- se_vec[keep]
    nms    <- nms[keep]
    term   <- term[keep]
  }

  # Wald intervals on the internal (estimation) scale.
  z         <- stats::qnorm((1 + level) / 2)
  estimates <- unname(coefs)
  conf_low  <- unname(coefs - z * se_vec)
  conf_high <- unname(coefs + z * se_vec)

  # Back-transform beta_k (k) AND log_gamma (gamma) rows to natural scale when
  # requested. log_sigma_u / beta0 stay on their internal scales (R6).
  if (report_space == "natural") {
    exp_pos <- which(nms %in% c("beta_k", "log_gamma"))
    if (length(exp_pos) > 0L) {
      estimates[exp_pos] <- exp(estimates[exp_pos])
      conf_low[exp_pos]  <- exp(conf_low[exp_pos])
      conf_high[exp_pos] <- exp(conf_high[exp_pos])
    }
  }

  tibble::tibble(
    term      = term,
    estimate  = estimates,
    conf.low  = conf_low,
    conf.high = conf_high,
    level     = level
  )
}


# --- summary / print ---

#' Summarize a structural choice discounting fit
#'
#' @param object A `beezdiscounting_choice` object.
#' @param report_space Scale for the fixed-effect (`beta_k`) estimates and
#'   standard errors in the coefficient table: `"natural"` (default; `k` via
#'   `exp()`), `"log10"`, `"internal"` or `"log"`. `statistic` and `p.value`
#'   are always computed on the estimation scale regardless of `report_space`.
#' @param ... Unused.
#' @return An object of class `"summary.beezdiscounting_choice"` with components:
#'   `call`, `model_class`, `backend`, `mode`, `equation`, `coefficients`,
#'   `variance_components`, `n_subjects`, `nobs`, `converged`, `logLik`,
#'   `AIC`, `BIC`, `notes`.
#'
#' @export
summary.beezdiscounting_choice <- function(object,
                                           report_space = c("natural", "log10", "internal", "log"),
                                           ...) {
  report_space <- match.arg(report_space)

  # Reuse tidy() for the fixed/shape coefficient rows (k + gamma + beta0),
  # applying the explicit gamma/identity-beta0 R6 handling consistently.
  td <- tidy(object, effects = "fixed", report_space = report_space)
  coefficients <- td[td$component %in% c("fixed", "shape"), , drop = FALSE]

  vc <- .dd_choice_variance_components(object)

  notes <- character(0)
  if (!object$converged) {
    notes <- c(notes, "WARNING: Model did not converge.")
  }
  if (isFALSE(object$se_available)) {
    notes <- c(notes, "Standard errors unavailable (sdreport failed); CIs will be NA.")
  }
  if (isFALSE(object$hessian_pd)) {
    notes <- c(notes,
      "Warning: Hessian not positive definite - standard errors may be unreliable.")
  }
  if (length(object$opt_warnings %||% character(0)) > 0L) {
    notes <- c(notes, sprintf(
      "Optimizer produced %d warning(s) during fitting.",
      length(object$opt_warnings)
    ))
  }
  if (!is.null(object$param_info$factors) && length(object$param_info$factors) > 0L) {
    notes <- c(notes,
      "Population k reflects the reference level. Use get_dd_param_emms() for per-group estimates.")
  }
  notes <- c(notes,
    "gamma is the choice-sensitivity (slope) parameter on the logit scale.")

  structure(
    list(
      call               = object$call,
      model_class        = "beezdiscounting_choice",
      backend            = "TMB_choice",
      mode               = object$param_info$mode,
      equation           = object$param_info$equation,
      coefficients       = coefficients,
      variance_components = vc,
      n_subjects         = object$param_info$n_subjects,
      nobs               = object$param_info$n_obs,
      converged          = object$converged,
      logLik             = object$loglik,
      AIC                = object$AIC,
      BIC                = object$BIC,
      notes              = notes
    ),
    class = "summary.beezdiscounting_choice"
  )
}


#' Print a structural choice discounting fit
#'
#' @param x A `beezdiscounting_choice` object.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#'
#' @export
print.beezdiscounting_choice <- function(x, ...) {
  cat("\nStructural SS-vs-LL Choice Discounting Model\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("Mode:", x$param_info$mode, "\n")
  cat("Equation:", x$param_info$equation, "\n")
  cat("Convergence:", ifelse(x$converged, "Yes", "No"), "\n")
  cat("Number of subjects:", x$param_info$n_subjects, "\n")
  cat("Number of observations:", x$param_info$n_obs, "\n")
  cat("Random effects:", x$param_info$n_random_effects, "(k ~ 1)\n")
  cat("Log-likelihood:", round(x$loglik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")

  cat("\nFixed Effects (log k):\n")
  co <- x$model$coefficients
  tn <- .dd_choice_term_names(x)
  fe <- co[tn$k_idx]
  names(fe) <- tn$term[tn$k_idx]
  print(round(fe, 4))

  cat("\nUse summary() for full results.\n")
  invisible(x)
}


#' Print a structural choice discounting model summary
#'
#' @param x A `summary.beezdiscounting_choice` object.
#' @param digits Number of significant digits for rounding.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#'
#' @export
print.summary.beezdiscounting_choice <- function(x, digits = 4, ...) {
  cat("\nStructural Choice Discounting Model Summary\n")
  cat(strrep("=", 50), "\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    cat(paste(deparse(x$call), collapse = "\n"), "\n\n")
  }
  cat("Mode:", x$mode, "\n")
  cat("Equation:", x$equation, "\n")
  cat("Backend:", x$backend, "\n")
  cat("Convergence:", ifelse(x$converged, "Yes", "No"), "\n")
  cat("Subjects:", x$n_subjects, " Observations:", x$nobs, "\n\n")

  scale_lbl <- unique(stats::na.omit(x$coefficients$estimate_scale[
    x$coefficients$component == "fixed"]))
  k_space <- switch(if (length(scale_lbl) == 1L) scale_lbl else "log",
    natural = "k",
    log10   = "log10 k",
    log     = "log k",
    "log k"
  )
  cat(sprintf("--- Fixed Effects (%s) ---\n", k_space))
  cd <- as.data.frame(x$coefficients[, c("term", "estimate", "std.error",
                                         "statistic", "p.value")])
  cd$estimate  <- round(cd$estimate, digits)
  cd$std.error <- round(cd$std.error, digits)
  cd$statistic <- round(cd$statistic, digits)
  cd$p.value   <- format.pval(cd$p.value, digits = 3)
  print(cd, row.names = FALSE)

  cat("\n--- Variance Components ---\n")
  vc <- x$variance_components
  vc$Estimate <- round(vc$Estimate, digits)
  print(vc, row.names = FALSE)

  cat("\n--- Fit Statistics ---\n")
  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "  BIC:", round(x$BIC, 2), "\n")

  if (length(x$notes) > 0L) {
    cat("\nNotes:\n")
    for (note in x$notes) cat("  *", note, "\n")
  }
  invisible(x)
}
