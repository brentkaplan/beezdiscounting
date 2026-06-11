# brms tier: S3 methods ---------------------------------------------------------
#
# Methods for beezdiscounting_brms, mirroring the beezdiscounting_tmb
# contracts (the exact 8-column tidy table, glance columns, confint shape)
# with draws-based uncertainty: report-space conversions transform the
# posterior draws and then summarize -- no delta method, no Wald statistics.

#' Fixed-effect (and shape) coefficient table from posterior draws
#' @noRd
.dd_brms_coef_table <- function(object, report_space = "natural") {
  report_space <- match.arg(report_space, c("natural", "log10", "internal", "log"))
  to <- if (report_space == "internal") "log" else report_space

  draws <- .dd_brms_draws_matrix(object)
  k_vars <- grep("^b_logk_", colnames(draws), value = TRUE)
  terms <- paste0("k:", colnames(object$formula_details$X))

  rows <- lapply(seq_along(k_vars), function(j) {
    dj <- as.numeric(draws[, k_vars[j]])
    tdraws <- .dd_brms_transform_draws(dj, to = to)
    tibble::tibble(
      term = terms[j],
      estimate = stats::median(tdraws),
      std.error = stats::sd(tdraws),
      statistic = NA_real_,
      p.value = NA_real_,
      component = "fixed",
      estimate_scale = to,
      term_display = .dd_term_display_space(terms[j], to)
    )
  })

  if (isTRUE(object$param_info$has_s)) {
    dj <- as.numeric(draws[, "b_logs_Intercept"])
    tdraws <- .dd_brms_transform_draws(dj, to = to)
    rows <- c(rows, list(tibble::tibble(
      term = "s",
      estimate = stats::median(tdraws),
      std.error = stats::sd(tdraws),
      statistic = NA_real_,
      p.value = NA_real_,
      component = "shape",
      estimate_scale = to,
      term_display = .dd_term_display_space("s", to)
    )))
  }

  dplyr::bind_rows(rows)
}

# --- basics ------------------------------------------------------------------------

#' @export
coef.beezdiscounting_brms <- function(object, ...) {
  object$model$coefficients
}

#' @export
fixef.beezdiscounting_brms <- function(object, ...) {
  object$model$coefficients
}

#' @export
ranef.beezdiscounting_brms <- function(object, ...) {
  draws <- .dd_brms_draws_matrix(object)
  ids <- object$param_info$subject_levels
  wanted <- paste0("r_id__logk[", ids, ",Intercept]")
  u <- apply(draws[, wanted, drop = FALSE], 2, stats::median)
  out <- data.frame(id = ids, r_logk = unname(u), stringsAsFactors = FALSE)
  sp <- object$subject_pars
  if (!is.null(sp)) {
    out$k <- sp$k[match(out$id, sp$id)]
  }
  out
}

#' @export
nobs.beezdiscounting_brms <- function(object, ...) {
  object$param_info$n_obs
}

#' @export
logLik.beezdiscounting_brms <- function(object, ...) {
  stop(
    "logLik() is not defined for Bayesian fits in this tier. ",
    "Use the stored $loo (or brms::loo(fit$brmsfit)) for model comparison.",
    call. = FALSE
  )
}

# --- broom --------------------------------------------------------------------------

#' Tidy a beezdiscounting_brms model
#'
#' The exact 8-column dd coefficient contract (`term`, `estimate`,
#' `std.error`, `statistic`, `p.value`, `component`, `estimate_scale`,
#' `term_display`). Estimates are posterior medians and `std.error`
#' posterior SDs of the report-space-transformed draws (exact; no delta
#' method); `statistic`/`p.value` are `NA` -- use `confint()` for interval
#' summaries. The shape row (`s`, two-parameter equations) carries
#' `component = "shape"`; variance rows (`effects = "ran_pars"`) mirror the
#' TMB reporting convention (log10-scale k RE SD; natural phi/sigma).
#'
#' @param x A `beezdiscounting_brms` object.
#' @param effects `"fixed"`, `"ran_pars"`, or both (default).
#' @param report_space `"natural"` (default), `"log10"`, `"internal"`, or
#'   `"log"`.
#' @param ... Unused.
#' @return A tibble.
#' @export
tidy.beezdiscounting_brms <- function(
  x,
  effects = c("fixed", "ran_pars"),
  report_space = c("natural", "log10", "internal", "log"),
  ...
) {
  effects <- match.arg(effects, several.ok = TRUE)
  report_space <- match.arg(report_space)

  result <- tibble::tibble(
    term = character(),
    estimate = numeric(),
    std.error = numeric(),
    statistic = numeric(),
    p.value = numeric(),
    component = character(),
    estimate_scale = character(),
    term_display = character()
  )

  if ("fixed" %in% effects) {
    result <- dplyr::bind_rows(result, .dd_brms_coef_table(x, report_space))
  }

  if ("ran_pars" %in% effects) {
    vc <- x$model$variance_components
    ran <- tibble::tibble(
      term = vc$Component,
      estimate = vc$Estimate,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      component = "variance",
      estimate_scale = vc$Scale,
      term_display = vc$Component
    )
    result <- dplyr::bind_rows(result, ran)
  }

  result
}

#' Glance at a beezdiscounting_brms model
#'
#' The dd glance columns with `backend = "brms"` and
#' `logLik`/`AIC`/`BIC = NA_real_` by design, plus the Bayesian comparison
#' and diagnostic currency (`elpd_loo`, `p_loo`, `looic`, `rhat_max`,
#' `ess_bulk_min`, `num_divergent`).
#'
#' @param x A `beezdiscounting_brms` object.
#' @param ... Unused.
#' @return A one-row tibble.
#' @export
glance.beezdiscounting_brms <- function(x, ...) {
  loo_est <- function(what) {
    if (is.null(x$loo)) {
      return(NA_real_)
    }
    est <- x$loo$estimates
    if (!what %in% rownames(est)) {
      return(NA_real_)
    }
    unname(est[what, "Estimate"])
  }

  tibble::tibble(
    model_class = "beezdiscounting_brms",
    backend = "brms",
    equation = x$param_info$equation,
    family = x$param_info$family,
    nobs = x$param_info$n_obs,
    n_subjects = x$param_info$n_subjects,
    n_random_effects = x$param_info$n_random_effects,
    converged = x$converged,
    logLik = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    elpd_loo = loo_est("elpd_loo"),
    p_loo = loo_est("p_loo"),
    looic = loo_est("looic"),
    rhat_max = x$mcmc_info$rhat_max,
    ess_bulk_min = x$mcmc_info$ess_bulk_min,
    num_divergent = x$mcmc_info$num_divergent
  )
}

#' @export
augment.beezdiscounting_brms <- function(x, ...) {
  fv <- fitted(x)
  out <- tibble::as_tibble(x$data)
  out$.fitted <- fv
  out$.resid <- x$data$y - fv
  disp <- .dd_brms_dispersion_sd(x)
  out$.std_resid <- out$.resid / disp
  out
}

#' Posterior-median residual scale: sigma (gaussian) or the beta-implied SD
#' at the fitted mean
#' @noRd
.dd_brms_dispersion_sd <- function(object) {
  draws <- .dd_brms_draws_matrix(object)
  if (object$param_info$family == "gaussian") {
    stats::median(as.numeric(draws[, "sigma"]))
  } else {
    phi <- stats::median(as.numeric(draws[, "phi"]))
    mu <- fitted(object)
    sqrt(mu * (1 - mu) / (1 + phi))
  }
}

# --- confint ------------------------------------------------------------------------

#' Credible intervals for a beezdiscounting_brms model
#'
#' Equal-tailed quantile credible intervals on the report-space-transformed
#' draws. Columns match `confint.beezdiscounting_tmb()`: `term`,
#' `estimate`, `conf.low`, `conf.high`, `level`.
#'
#' @param object A `beezdiscounting_brms` object.
#' @param parm Optional terms (display, e.g. `"k:(Intercept)"`/`"s"`, or
#'   TMB coefficient names `"beta_k"`/`"log_s"`).
#' @param level Credible level (default 0.95).
#' @param report_space Reporting scale.
#' @param ... Unused.
#' @return A tibble.
#' @export
confint.beezdiscounting_brms <- function(
  object,
  parm = NULL,
  level = 0.95,
  report_space = c("natural", "log10", "internal", "log"),
  ...
) {
  report_space <- match.arg(report_space)
  to <- if (report_space == "internal") "log" else report_space
  alpha2 <- (1 - level) / 2

  draws <- .dd_brms_draws_matrix(object)
  k_vars <- grep("^b_logk_", colnames(draws), value = TRUE)
  terms <- paste0("k:", colnames(object$formula_details$X))
  tmb_names <- rep("beta_k", length(k_vars))
  vars <- k_vars
  if (isTRUE(object$param_info$has_s)) {
    vars <- c(vars, "b_logs_Intercept")
    terms <- c(terms, "s")
    tmb_names <- c(tmb_names, "log_s")
  }

  rows <- lapply(seq_along(vars), function(j) {
    tdraws <- .dd_brms_transform_draws(as.numeric(draws[, vars[j]]), to = to)
    tibble::tibble(
      term = terms[j],
      estimate = stats::median(tdraws),
      conf.low = unname(stats::quantile(tdraws, alpha2)),
      conf.high = unname(stats::quantile(tdraws, 1 - alpha2)),
      level = level
    )
  })
  out <- dplyr::bind_rows(rows)

  if (!is.null(parm)) {
    display <- vapply(terms, .dd_term_display_space, character(1), report_space = to)
    keep <- out$term %in% parm | tmb_names %in% parm | display %in% parm
    out <- out[keep, , drop = FALSE]
  }
  out
}

# --- predict / fitted / residuals ----------------------------------------------------

#' Predict from a beezdiscounting_brms model
#'
#' `type = "response"` summarizes `brms::posterior_epred()` draws (the
#' indifference-proportion mean) to the posterior median with equal-tailed
#' intervals; `level = "subject"` conditions on the subject random effects,
#' `"population"` sets them to zero. `type = "parameters"` returns the
#' per-subject k summaries.
#'
#' @param object A `beezdiscounting_brms` object.
#' @param newdata Optional data frame (defaults to the model frame). The
#'   Rachlin guard columns (`xzero`/`xsafe`) are derived automatically.
#' @param type `"response"` or `"parameters"`.
#' @param level `"subject"` (default) or `"population"`.
#' @param probs Interval probabilities.
#' @param draws If `TRUE`, return the epred draws matrix.
#' @param ... Unused.
#' @return A tibble (or draws matrix).
#' @export
predict.beezdiscounting_brms <- function(
  object,
  newdata = NULL,
  type = c("response", "parameters"),
  level = "subject",
  probs = c(0.025, 0.975),
  draws = FALSE,
  ...
) {
  type <- match.arg(type)
  level <- match.arg(level, c("subject", "population"))

  if (type == "parameters") {
    return(object$subject_pars)
  }

  if (!is.null(newdata) && object$param_info$equation == "rachlin" &&
    !all(c("xzero", "xsafe") %in% names(newdata))) {
    newdata$xzero <- as.numeric(newdata$x == 0)
    newdata$xsafe <- ifelse(newdata$x == 0, 1, newdata$x)
  }

  re_formula <- if (level == "population") NA else NULL
  ep <- brms::posterior_epred(
    object$brmsfit,
    newdata = newdata, re_formula = re_formula
  )
  if (draws) {
    return(ep)
  }

  base <- tibble::as_tibble(if (is.null(newdata)) object$data else newdata)
  base$.fitted <- apply(ep, 2, stats::median)
  base$.lower <- apply(ep, 2, stats::quantile, probs = probs[1])
  base$.upper <- apply(ep, 2, stats::quantile, probs = probs[2])
  base
}

#' @export
fitted.beezdiscounting_brms <- function(object, ...) {
  ep <- brms::posterior_epred(object$brmsfit, re_formula = NULL)
  apply(ep, 2, stats::median)
}

#' @export
residuals.beezdiscounting_brms <- function(
  object,
  type = c("response", "pearson"),
  ...
) {
  type <- match.arg(type)
  res <- object$data$y - fitted(object)
  if (type == "pearson") {
    res <- res / .dd_brms_dispersion_sd(object)
  }
  res
}

# --- print / summary ----------------------------------------------------------------

#' @export
print.beezdiscounting_brms <- function(x, digits = 4, ...) {
  cat("Bayesian Mixed-Effects Discounting Model (brms)\n")
  cat(strrep("=", 50), "\n")
  cat("Equation:", x$param_info$equation, " Family:", x$param_info$family, "\n")
  cat(
    "Chains:", x$mcmc_info$chains,
    " Iterations:", x$mcmc_info$iter,
    " (warmup", paste0(x$mcmc_info$warmup, ")"), "\n"
  )
  cat(
    "Subjects:", x$param_info$n_subjects,
    " Observations:", x$param_info$n_obs, "\n"
  )
  cat("Converged:", ifelse(isTRUE(x$converged), "Yes", "No"), "\n")
  cat("\nCoefficients (posterior medians, natural scale):\n")
  ct <- .dd_brms_coef_table(x, "natural")
  print(
    as.data.frame(ct[, c("term_display", "estimate", "std.error")]),
    digits = digits, row.names = FALSE
  )
  cat("\nUse summary(), tidy(), or confint() for full posterior summaries.\n")
  invisible(x)
}

#' Summarize a beezdiscounting_brms model
#'
#' @param object A `beezdiscounting_brms` object.
#' @param report_space Reporting scale for the coefficient table.
#' @param ... Unused.
#' @return A `summary.beezdiscounting_brms` list.
#' @export
summary.beezdiscounting_brms <- function(
  object,
  report_space = c("natural", "log10", "internal", "log"),
  ...
) {
  report_space <- match.arg(report_space)
  structure(
    list(
      equation = object$param_info$equation,
      family = object$param_info$family,
      backend = "brms",
      converged = object$converged,
      n_subjects = object$param_info$n_subjects,
      nobs = object$param_info$n_obs,
      coefficients = .dd_brms_coef_table(object, report_space),
      variance_components = object$model$variance_components,
      mcmc_info = object$mcmc_info,
      autoscale_info = object$autoscale_info,
      report_space = report_space
    ),
    class = "summary.beezdiscounting_brms"
  )
}

#' @export
print.summary.beezdiscounting_brms <- function(x, digits = 4, ...) {
  cat("\nBayesian Mixed-Effects Discounting Model Summary (brms)\n")
  cat(strrep("=", 50), "\n\n")
  cat("Equation:", x$equation, " Family:", x$family, "\n")
  cat("Backend:", x$backend, "\n")
  cat("Converged:", ifelse(isTRUE(x$converged), "Yes", "No"), "\n")
  cat("Subjects:", x$n_subjects, " Observations:", x$nobs, "\n")
  cat(sprintf(
    "Diagnostics: rhat_max = %.3f, divergences = %d, min bulk ESS = %d\n\n",
    x$mcmc_info$rhat_max,
    x$mcmc_info$num_divergent,
    round(x$mcmc_info$ess_bulk_min)
  ))
  cat("Coefficients (posterior median / SD,", x$report_space, "scale):\n")
  print(
    as.data.frame(
      x$coefficients[, c("term_display", "estimate", "std.error")]
    ),
    digits = digits, row.names = FALSE
  )
  cat("\nVariance components:\n")
  print(x$variance_components, digits = digits, row.names = FALSE)
  invisible(x)
}
