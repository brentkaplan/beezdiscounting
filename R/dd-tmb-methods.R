# ==============================================================================
# S3 Methods for beezdiscounting_tmb Objects
# ==============================================================================

#' Build display term names from a beezdiscounting_tmb coefficient vector
#'
#' Maps raw optimizer names (`beta_k`, `log_sigma_u`, `log_phi` / `log_sigma_e`)
#' to readable display names. `beta_k` columns become `k:<design colname>`
#' using `formula_details$X` colnames; every other coefficient keeps its raw
#' name. Shared by `tidy()`, `summary()`, `confint()`.
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param nms Character vector of raw parameter names (default
#'   `names(object$model$coefficients)`).
#' @return List with `term` (display names), `k_idx` (beta_k positions),
#'   `other_idx` (non-beta positions).
#' @keywords internal
.dd_tmb_build_term_names <- function(object, nms = NULL) {
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

  list(
    term = term,
    k_idx = k_idx,
    other_idx = as.integer(other_idx)
  )
}


#' Named standard-error vector aligned to the coefficient vector
#'
#' Pulls fixed-effect SEs from the sdreport. The optimizer parameter vector
#' (`beta_k`, `log_sigma_u`, `log_aux`) is what `sdreport$par.fixed` /
#' `sdreport$cov.fixed` cover; entries are renamed to match
#' `names(object$model$coefficients)` (i.e. `log_aux` already displayed as
#' `log_phi` / `log_sigma_e` in the fit object). Returns `NA` SEs when the
#' sdreport is unavailable.
#'
#' @param object A `beezdiscounting_tmb` object.
#' @return Named numeric vector parallel to `object$model$coefficients`.
#' @keywords internal
.dd_tmb_model_se <- function(object) {
  co <- object$model$coefficients
  if (!is.null(object$model$se)) {
    se <- object$model$se
    # Defensive: align to coefficient order/names.
    return(stats::setNames(unname(se)[match(names(co), names(se))], names(co)))
  }
  sdr <- object$sdr
  if (is.null(sdr) || isFALSE(object$se_available)) {
    return(stats::setNames(rep(NA_real_, length(co)), names(co)))
  }
  sd_fixed <- sqrt(diag(as.matrix(sdr$cov.fixed)))
  raw_nms <- names(sdr$par.fixed)
  # `log_aux` is the optimizer name; the fit object renames it to log_phi /
  # log_sigma_e. Translate so the SE vector lines up with coefficients.
  aux_name <- intersect(c("log_phi", "log_sigma_e"), names(co))
  raw_nms[raw_nms == "log_aux"] <- if (length(aux_name)) aux_name[1] else "log_aux"
  se <- stats::setNames(rep(NA_real_, length(co)), names(co))
  hit <- match(names(co), raw_nms)
  se[!is.na(hit)] <- sd_fixed[hit[!is.na(hit)]]
  se
}


# --- logLik / AIC / BIC / nobs ---

#' @export
logLik.beezdiscounting_tmb <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- length(object$opt$par)
  attr(ll, "nobs") <- object$param_info$n_obs
  class(ll) <- "logLik"
  ll
}

#' @export
AIC.beezdiscounting_tmb <- function(object, ..., k = 2) {
  if (k != 2) {
    nll <- -object$loglik
    n_params <- length(object$opt$par)
    return(2 * nll + k * n_params)
  }
  object$AIC
}

#' @export
BIC.beezdiscounting_tmb <- function(object, ...) {
  object$BIC
}

#' @export
nobs.beezdiscounting_tmb <- function(object, ...) {
  object$param_info$n_obs %||% nrow(object$data)
}


# --- coef / fixef / ranef ---

#' Extract coefficients from a TMB discounting model
#'
#' Returns the optimizer's flat named parameter vector: `beta_k` (one entry per
#' fixed-effect design column, on the log-k scale), `log_sigma_u`, and the
#' auxiliary parameter (`log_phi` for `family = "sltb"`, `log_sigma_e` for
#' `family = "gaussian"`). This is the numeric escape hatch consumed by tooling
#' such as `car::deltaMethod`.
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return Named numeric vector.
#' @export
coef.beezdiscounting_tmb <- function(object, ...) {
  object$model$coefficients
}

#' Extract fixed effects from a TMB discounting model
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return Named numeric vector (identical to `coef()`).
#' @importFrom nlme fixef
#' @export
fixef.beezdiscounting_tmb <- function(object, ...) {
  coef(object)
}

#' Extract subject-level random effects from a TMB discounting model
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return Data frame with `id`, the standardized random-intercept deviate
#'   `u_i` (such that `log k_i = X beta + sigma_u * u_i`), and the resolved
#'   per-subject discount rate `k`. There is **no** `phi` column: phi is
#'   population-level in the MVP, not a subject-level parameter.
#' @importFrom nlme ranef
#' @export
ranef.beezdiscounting_tmb <- function(object, ...) {
  sp <- object$subject_pars
  keep <- intersect(c("id", "u_i", "k"), names(sp))
  out <- sp[, keep, drop = FALSE]
  rownames(out) <- NULL
  out
}


# --- predict ---

#' Discounting function value mu = E[y] for a vector of k and x
#'
#' Mazur: `mu = 1 / (1 + k * x)`; exponential: `mu = exp(-k * x)`. Guards mu
#' to `[1e-6, 1 - 1e-6]` to match the C++ template bounds.
#'
#' @param k Numeric vector of discount rates.
#' @param x Numeric vector of delays (same length as `k`).
#' @param equation Character, `"mazur"` or `"exponential"`.
#' @return Numeric vector of mu values clamped to `[1e-6, 1-1e-6]`.
#' @keywords internal
.dd_discount_mu <- function(k, x, equation) {
  mu <- switch(equation,
    mazur = 1 / (1 + k * x),
    exponential = exp(-k * x),
    stop("unknown equation '", equation, "'", call. = FALSE)
  )
  pmin(pmax(mu, 1e-6), 1 - 1e-6)
}


#' Rebuild the per-row log-k linear predictor from newdata
#'
#' Reconstructs the fixed-effect design from `newdata` using the stored RHS
#' formula and contrasts, multiplies by `beta_k`, and (for
#' `level = "subject"`) adds each subject's `sigma_u * u_i` deviate looked up
#' by id. Returns per-row `k` on the natural scale.
#'
#' Unseen factor levels in `newdata` are rejected with a clear error rather
#' than silently zero-padded (R2 fix).
#'
#' @param object A `beezdiscounting_tmb` fit.
#' @param newdata Data frame with the x var, factor/covariate cols, and (for
#'   `level = "subject"`) the id var.
#' @param level `"subject"` or `"population"`.
#' @return Numeric vector of per-row `k`, length `nrow(newdata)`.
#' @keywords internal
.dd_tmb_predict_k <- function(object, newdata, level = c("subject", "population")) {
  level <- match.arg(level)
  pinfo  <- object$param_info
  coefs  <- object$model$coefficients
  beta_k <- unname(coefs[names(coefs) == "beta_k"])

  # Rebuild the newdata design through the SAME route as the fit: the stored
  # one-sided rhs formula + stored per-factor contrasts (R2 fix).  Passing
  # contrasts.arg pins the contrast coding so the columns match the fit; we do
  # NOT silently zero-pad missing columns (that would mask a level mismatch).
  rhs           <- object$formula_details$rhs %||% stats::as.formula("~ 1")
  contrasts_arg <- object$formula_details$contrasts

  # ERROR (do not zero-pad) on unseen factor levels.
  for (f in pinfo$factors) {
    if (!is.null(f) && nzchar(f) && f %in% names(newdata) &&
          f %in% names(object$data)) {
      fit_levels <- levels(as.factor(object$data[[f]]))
      nd_levels  <- unique(as.character(newdata[[f]]))
      unseen     <- setdiff(nd_levels, fit_levels)
      if (length(unseen) > 0L) {
        cli::cli_abort(c(
          "newdata factor {.field {f}} has level(s) not seen in the fit: {.val {unseen}}.",
          i = "predict() cannot encode unseen factor levels."
        ))
      }
      # Re-level to the fitted factor so model.matrix builds the SAME columns.
      newdata[[f]] <- factor(newdata[[f]], levels = fit_levels)
    }
  }

  Xnew   <- stats::model.matrix(rhs, data = newdata,
                                contrasts.arg = contrasts_arg)
  fit_cn <- colnames(object$formula_details$X)
  if (!identical(colnames(Xnew), fit_cn)) {
    if (!all(fit_cn %in% colnames(Xnew))) {
      cli::cli_abort(c(
        "newdata design columns do not match the fitted design.",
        i = "Expected: {.val {fit_cn}}; got: {.val {colnames(Xnew)}}."
      ))
    }
    Xnew <- Xnew[, fit_cn, drop = FALSE]
  }

  eta <- as.numeric(Xnew %*% beta_k)

  if (level == "subject") {
    id_var <- pinfo$id_var
    if (is.null(newdata[[id_var]])) {
      cli::cli_abort(c(
        "{.code level = \"subject\"} needs the id column {.val {id_var}} in {.arg newdata}.",
        i = "Use {.code level = \"population\"} for the random-effects-at-zero curve."
      ))
    }
    sp      <- object$subject_pars
    sigma_u <- exp(unname(coefs[["log_sigma_u"]]))
    u_by_id <- stats::setNames(sp$u_i, as.character(sp$id))
    u_row   <- u_by_id[as.character(newdata[[id_var]])]
    if (anyNA(u_row)) {
      unknown_ids <- unique(as.character(newdata[[id_var]])[is.na(u_row)])
      cli::cli_abort(c(
        "newdata contains id{?s} not present in the fit: {.val {unknown_ids}}.",
        i = "Use {.code level = \"population\"} for out-of-sample predictions."
      ))
    }
    eta <- eta + sigma_u * unname(u_row)
  }

  exp(eta)
}


#' Predict from a TMB mixed-effects discounting model
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param newdata Optional data frame. `NULL` uses the fitting data.
#' @param type `"response"` (fitted indifference proportions on `(0,1)`) or
#'   `"parameters"` (the per-subject parameter tibble).
#' @param level For `type = "response"`: `"subject"` (default; conditions on
#'   each subject's estimated random intercept, requires the id column) and/or
#'   `"population"` (random effects set to zero — the population-mean curve;
#'   no id column needed). Pass `c("population", "subject")` for both columns
#'   side-by-side. A numeric nlme-style level is rejected with an error.
#' @param ... Unused.
#'
#' @return
#' - `type = "parameters"`: the subject-parameter tibble (one row per subject).
#' - `type = "response"`, `level = "subject"` only: `newdata` as a tibble plus
#'   a `.fitted` column.
#' - `type = "response"`, `"population"` requested: `newdata` plus
#'   `predict.fixed` (and `predict.id` when `"subject"` is also requested),
#'   matching the `nlme::predict.lme(level = 0:1)` column-name convention so
#'   nlme-based plotting code runs unchanged.
#'
#' @examples
#' \donttest{
#' dd <- .simulate_dd_ip_mixed(n_subjects = 20, seed = 1)
#' fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb",
#'                   random_effects = k ~ 1, verbose = 0)
#'
#' # Subject-conditional fitted values
#' head(predict(fit, type = "response"))
#'
#' # Population-mean curve at specific delays (no id column needed)
#' nd <- data.frame(x = c(7, 30, 180, 365))
#' predict(fit, newdata = nd, level = "population")
#'
#' # Per-subject parameters
#' head(predict(fit, type = "parameters"))
#' }
#'
#' @export
predict.beezdiscounting_tmb <- function(object,
                                        newdata  = NULL,
                                        type     = c("response", "parameters"),
                                        level    = "subject",
                                        ...) {
  type <- match.arg(type)
  # Reject numeric nlme-style level (0/1) with a clear error.  match.arg()
  # accepts only character; a numeric level produces an opaque error, so we
  # check first.
  if (!is.character(level)) {
    cli::cli_abort(
      c("{.arg level} must be a character vector.",
        i = "{.arg level} should be one of {.val \"subject\"} or {.val \"population\"}."),
      call = NULL
    )
  }
  level <- match.arg(level, c("subject", "population"), several.ok = TRUE)

  if (type == "parameters") {
    return(tibble::as_tibble(object$subject_pars))
  }

  # type == "response"
  equation <- object$param_info$equation
  x_var    <- object$param_info$x_var
  if (is.null(newdata)) newdata <- object$data
  out <- tibble::as_tibble(newdata)
  x   <- newdata[[x_var]]

  # Single "subject" level: historical .fitted column name for backward compat.
  if (identical(level, "subject")) {
    k_row       <- .dd_tmb_predict_k(object, newdata, level = "subject")
    out$.fitted <- .dd_discount_mu(k_row, x, equation)
    return(out)
  }

  # "population" and/or both: use nlme-style column names.
  if ("population" %in% level) {
    k_pop            <- .dd_tmb_predict_k(object, newdata, level = "population")
    out$predict.fixed <- .dd_discount_mu(k_pop, x, equation)
  }
  if ("subject" %in% level) {
    k_sub          <- .dd_tmb_predict_k(object, newdata, level = "subject")
    out$predict.id <- .dd_discount_mu(k_sub, x, equation)
  }
  out
}


# --- fitted / residuals / augment ---

#' Per-row response SD on the [0,1] scale for standardized residuals
#'
#' Gaussian: constant `sigma_e = exp(log_sigma_e)`. SLT-beta: the
#' delta-method SLT SD `s * sqrt(mu * (1 - mu) / (phi + 1))` (the SLT
#' variance at `s ≈ 1`), so residuals near the bounds are down-weighted.
#'
#' @param object A `beezdiscounting_tmb` fit.
#' @param mu Numeric vector of fitted mu values (from `.dd_discount_mu`).
#' @return Numeric vector of per-row response SDs, same length as `mu`.
#' @keywords internal
.dd_tmb_response_sd <- function(object, mu) {
  coefs  <- object$model$coefficients
  family <- object$param_info$family
  if (family == "gaussian") {
    return(rep(exp(coefs[["log_sigma_e"]]), length(mu)))
  }
  # sltb: delta-method SD with s ≈ 1 (s is the SLT scale parameter; the MVP
  # fixes it implicitly at 1).  The SLT variance is s^2 * mu*(1-mu)/(phi+1),
  # so SD = s * sqrt(mu*(1-mu)/(phi+1)).
  s   <- 1
  phi <- exp(coefs[["log_phi"]])
  s * sqrt(mu * (1 - mu) / (phi + 1))
}


#' Shared fitted-values / residuals back end
#'
#' Calls `predict()` at the requested level and returns `.fitted`, `.resid`,
#' and the data frame used (either `object$data` or `newdata`).
#'
#' @param object A `beezdiscounting_tmb` fit.
#' @param newdata Optional data frame; `NULL` → training data.
#' @param level `"subject"` (default) or `"population"`.
#' @return List with `.fitted`, `.resid`, and `data`.
#' @keywords internal
.dd_tmb_fitted_resid <- function(object, newdata = NULL,
                                 level = c("subject", "population")) {
  level      <- match.arg(level)
  data_used  <- if (is.null(newdata)) object$data else newdata
  pred       <- predict(object, newdata = data_used,
                        type = "response", level = level)
  fitted_col <- if (level == "population") "predict.fixed" else ".fitted"
  fitted_vals <- pred[[fitted_col]]
  y_obs       <- data_used[[object$param_info$y_var]]
  list(.fitted = fitted_vals, .resid = y_obs - fitted_vals, data = data_used)
}


#' Fitted values for a beezdiscounting_tmb fit
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param level `"subject"` (default) or `"population"`.
#' @param ... Unused.
#' @return Numeric vector of fitted indifference proportions, length `nobs(object)`.
#'
#' @examples
#' \donttest{
#' dd <- .simulate_dd_ip_mixed(n_subjects = 20, seed = 1)
#' fit <- fit_dd_tmb(dd, verbose = 0)
#' head(fitted(fit))
#' }
#'
#' @export
fitted.beezdiscounting_tmb <- function(object,
                                       level = c("subject", "population"),
                                       ...) {
  level <- match.arg(level)
  .dd_tmb_fitted_resid(object, level = level)$.fitted
}


#' Residuals for a beezdiscounting_tmb fit
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param type `"response"` (default; raw `y - fitted`) or `"pearson"` (divided
#'   by the per-row response SD; equivalent to `.std_resid` in `augment()`).
#' @param level `"subject"` (default) or `"population"`.
#' @param ... Unused.
#' @return Numeric vector of residuals, length `nobs(object)`.
#'
#' @examples
#' \donttest{
#' dd <- .simulate_dd_ip_mixed(n_subjects = 20, seed = 1)
#' fit <- fit_dd_tmb(dd, verbose = 0)
#' head(residuals(fit))
#' head(residuals(fit, type = "pearson"))
#' }
#'
#' @export
residuals.beezdiscounting_tmb <- function(object,
                                          type  = c("response", "pearson"),
                                          level = c("subject", "population"),
                                          ...) {
  type  <- match.arg(type)
  level <- match.arg(level)
  fr    <- .dd_tmb_fitted_resid(object, level = level)
  if (type == "response") return(fr$.resid)
  fr$.resid / .dd_tmb_response_sd(object, fr$.fitted)
}


#' Augment a beezdiscounting_tmb model
#'
#' Returns the data used for fitting (or `newdata`) as a tibble with three
#' diagnostic columns appended:
#' - `.fitted`: subject-conditional fitted indifference proportion (clamped to
#'   `(0, 1)`).
#' - `.resid`: raw residual `y - .fitted` on the response scale.
#' - `.std_resid`: Pearson (standardized) residual — `.resid` divided by the
#'   per-row response SD.  For `family = "gaussian"` the SD is the constant
#'   `sigma_e`; for `family = "sltb"` it is the delta-method SLT SD
#'   `sqrt(mu * (1 - mu) / (phi + 1))`.
#'
#' @param x A `beezdiscounting_tmb` object.
#' @param newdata Optional data frame. `NULL` uses the fitting data.
#' @param ... Unused.
#' @return A tibble with the same rows as the data plus `.fitted`, `.resid`,
#'   and `.std_resid`.
#'
#' @examples
#' \donttest{
#' dd <- .simulate_dd_ip_mixed(n_subjects = 20, seed = 1)
#' fit <- fit_dd_tmb(dd, verbose = 0)
#' head(augment(fit))
#' }
#'
#' @importFrom generics augment
#' @export
augment.beezdiscounting_tmb <- function(x, newdata = NULL, ...) {
  fr          <- .dd_tmb_fitted_resid(x, newdata = newdata, level = "subject")
  out         <- tibble::as_tibble(fr$data)
  out$.fitted   <- fr$.fitted
  out$.resid    <- fr$.resid
  out$.std_resid <- fr$.resid / .dd_tmb_response_sd(x, fr$.fitted)
  out
}
