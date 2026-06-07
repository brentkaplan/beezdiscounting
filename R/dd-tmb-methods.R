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
