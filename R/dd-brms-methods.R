# brms tier: S3 methods ---------------------------------------------------------
#
# Methods for beezdiscounting_brms, mirroring the beezdiscounting_tmb
# contracts (the exact 8-column tidy table, glance columns, confint shape)
# with draws-based uncertainty: report-space conversions transform the
# posterior draws and then summarize -- no delta method, no Wald statistics.

#' Fixed-effect (and shape) coefficient table from posterior draws
#' @noRd
.dd_brms_coef_table <- function(object, report_space = "natural") {
  report_space <- match.arg(
    report_space,
    c("natural", "log10", "internal", "log")
  )
  to <- if (report_space == "internal") "log" else report_space

  draws <- .dd_brms_draws_matrix(object)
  k_vars <- .dd_brms_logk_draw_vars(object, colnames(draws))
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
    rows <- c(
      rows,
      list(tibble::tibble(
        term = "s",
        estimate = stats::median(tdraws),
        std.error = stats::sd(tdraws),
        statistic = NA_real_,
        p.value = NA_real_,
        component = "shape",
        estimate_scale = to,
        term_display = .dd_term_display_space("s", to)
      ))
    )
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
  if (isTRUE(object$param_info$n_random_effects == 2L)) {
    # k + phi ~ 1 also models a per-subject precision RE; report it alongside
    # the log k deviate (and the natural-scale phi from subject_pars).
    rphi <- paste0("r_id__phi[", ids, ",Intercept]")
    out$r_phi <- unname(apply(draws[, rphi, drop = FALSE], 2, stats::median))
    if (!is.null(sp) && !is.null(sp$phi)) {
      out$phi <- sp$phi[match(out$id, sp$id)]
    }
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
  # one epred pass shared by the fitted values and the dispersion scale
  # (Codex whole-branch review R1/R4)
  ep <- brms::posterior_epred(x$brmsfit, re_formula = NULL)
  fv <- apply(ep, 2, stats::median)
  out <- tibble::as_tibble(x$data)
  out$.fitted <- fv
  out$.resid <- x$data$y - fv
  out$.std_resid <- out$.resid / .dd_brms_dispersion_sd(x, epred_draws = ep)
  out
}

#' Posterior-median residual scale: sigma (gaussian) or the beta-implied SD
#'
#' Draws-first throughout (Codex whole-branch review R1): for the beta
#' family the per-observation SD draws sqrt(mu_d (1 - mu_d) / (1 + phi_d))
#' are computed per draw and THEN summarized -- never sqrt() of summarized
#' mu/phi (the transform is nonlinear, so the order matters).
#' @noRd
.dd_brms_dispersion_sd <- function(object, epred_draws = NULL) {
  draws <- .dd_brms_draws_matrix(object)
  if (object$param_info$family == "gaussian") {
    stats::median(as.numeric(draws[, "sigma"]))
  } else {
    if (is.null(epred_draws)) {
      epred_draws <- brms::posterior_epred(object$brmsfit, re_formula = NULL)
    }
    if (isTRUE(object$param_info$n_random_effects == 2L)) {
      # k + phi ~ 1: phi is predicted, so it varies per observation. Use the
      # per-draw, per-observation precision (log-link inverse, REs included),
      # which is matrix-aligned with the epred draws -- there is no scalar
      # `phi` draw for a 2-RE fit.
      phi_mat <- brms::posterior_linpred(
        object$brmsfit,
        dpar = "phi",
        transform = TRUE,
        re_formula = NULL
      )
      sd_draws <- sqrt(epred_draws * (1 - epred_draws) / (1 + phi_mat))
    } else {
      phi_draws <- as.numeric(draws[, "phi"])
      # phi recycles down the columns (draw-aligned: column length == n_draws)
      if (nrow(epred_draws) != length(phi_draws)) {
        stop(
          "Internal error: epred draw count (",
          nrow(epred_draws),
          ") does not match phi draw count (",
          length(phi_draws),
          ").",
          call. = FALSE
        )
      }
      sd_draws <- sqrt(epred_draws * (1 - epred_draws) / (1 + phi_draws))
    }
    apply(sd_draws, 2, stats::median)
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
  k_vars <- .dd_brms_logk_draw_vars(object, colnames(draws))
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
    display <- vapply(
      terms,
      .dd_term_display_space,
      character(1),
      report_space = to
    )
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

  if (
    !is.null(newdata) &&
      object$param_info$equation == "rachlin" &&
      !all(c("xzero", "xsafe") %in% names(newdata))
  ) {
    newdata$xzero <- as.numeric(newdata$x == 0)
    newdata$xsafe <- ifelse(newdata$x == 0, 1, newdata$x)
  }

  re_formula <- if (level == "population") NA else NULL
  ep <- brms::posterior_epred(
    object$brmsfit,
    newdata = newdata,
    re_formula = re_formula
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
  ep <- brms::posterior_epred(object$brmsfit, re_formula = NULL)
  res <- object$data$y - apply(ep, 2, stats::median)
  if (type == "pearson") {
    res <- res / .dd_brms_dispersion_sd(object, epred_draws = ep)
  }
  res
}

# --- EMMs and comparisons (draws-based; TICKET-041) -------------------------------

#' Draws-based EMMs of k for brms fits (delegated from get_dd_param_emms)
#'
#' Same reference grid as the TMB path (`.dd_build_emm_ref_grid()`); each
#' cell's linear predictor is computed per posterior draw and summarized as
#' the posterior median with equal-tailed quantile intervals (no delta
#' method). Columns mirror the TMB output: `level`, `k`, `k_log`,
#' `std.error` (log-scale posterior SD), `conf.low`, `conf.high`.
#' @noRd
.dd_brms_param_emms <- function(
  fit,
  factors_in_emm = NULL,
  at = NULL,
  ci_level = 0.95
) {
  .dd_validate_at(fit, at)
  grid <- .dd_build_emm_ref_grid(
    fit,
    at = at,
    factors_in_emm = factors_in_emm,
    validate = FALSE
  )
  draws <- .dd_brms_draws_matrix(fit)
  b <- as.matrix(
    draws[, .dd_brms_logk_draw_vars(fit, colnames(draws)), drop = FALSE]
  )
  alpha2 <- (1 - ci_level) / 2

  summarize_lin <- function(lin, label) {
    tibble::tibble(
      level = label,
      k = stats::median(exp(lin)),
      k_log = stats::median(lin),
      std.error = stats::sd(lin),
      conf.low = unname(stats::quantile(exp(lin), alpha2)),
      conf.high = unname(stats::quantile(exp(lin), 1 - alpha2))
    )
  }

  if (isTRUE(grid$is_intercept_only)) {
    return(summarize_lin(as.numeric(b[, 1]), "(Intercept)"))
  }

  ref_X <- grid$ref_X
  if (ncol(ref_X) != ncol(b)) {
    stop(
      "Reference-grid design and posterior coefficient draws disagree.",
      call. = FALSE
    )
  }
  cell_lin <- b %*% t(ref_X)
  dplyr::bind_rows(lapply(seq_len(nrow(ref_X)), function(i) {
    summarize_lin(
      as.numeric(cell_lin[, i]),
      .dd_brms_grid_label(grid, i, grid$use_factors)
    )
  }))
}

#' Reference-grid cell label (TMB label vocabulary)
#' @noRd
.dd_brms_grid_label <- function(grid, i, fs) {
  if (length(fs) > 0L) {
    paste(
      vapply(
        fs,
        function(f) {
          paste0(f, "=", grid$level_combos[[f]][i])
        },
        character(1)
      ),
      collapse = ", "
    )
  } else if (length(grid$cov_names) > 0L) {
    paste(
      vapply(
        grid$cov_names,
        function(cv) {
          paste0(cv, "=", grid$level_combos[[cv]][i])
        },
        character(1)
      ),
      collapse = ", "
    )
  } else {
    "(Intercept)"
  }
}

#' Posterior probability of direction, tie-aware (Codex 041-R1)
#' @noRd
.dd_brms_post_prob <- function(d) {
  max(mean(d > 0), mean(d < 0)) + 0.5 * mean(d == 0)
}

#' Draws-based k contrasts for brms fits (delegated from get_dd_comparisons)
#'
#' Per-draw differences of reference-grid linear predictors, summarized on
#' the log10 scale and as natural-scale ratios with quantile credible
#' intervals. No multiplicity adjustment (the joint posterior already
#' encodes contrast dependence): `statistic`/`df`/`p.value` are `NA` and
#' `post.prob` reports the posterior probability of direction. The return
#' mirrors the TMB `beezdiscounting_comparison` structure (backend
#' `"brms"`), so `print()`/`tidy()` work unchanged.
#' @noRd
.dd_brms_comparisons <- function(
  fit,
  compare_specs = NULL,
  contrast_type = c("pairwise", "trt.vs.ctrl"),
  contrast_by = NULL,
  adjust = "none",
  at = NULL,
  ci_level = 0.95,
  report_ratios = TRUE
) {
  contrast_type <- match.arg(contrast_type)
  if (!identical(adjust, "none")) {
    warning(
      "Multiplicity adjustment ('",
      adjust,
      "') does not apply to posterior ",
      "summaries; reporting unadjusted quantile intervals (see post.prob).",
      call. = FALSE
    )
  }

  fitted_factors <- fit$param_info$factors
  fitted_factors <- fitted_factors[
    !is.na(fitted_factors) & nzchar(fitted_factors)
  ]
  factors_in_emm <- NULL
  if (!is.null(compare_specs)) {
    if (!inherits(compare_specs, "formula")) {
      stop(
        "`compare_specs` must be a one-sided formula (e.g. ~ group).",
        call. = FALSE
      )
    }
    factors_in_emm <- all.vars(compare_specs)
    bad <- setdiff(factors_in_emm, fitted_factors)
    if (length(bad) > 0L) {
      stop(
        "`compare_specs` names factor(s) not in the fit: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (!is.null(contrast_by)) {
    bad_by <- setdiff(contrast_by, fitted_factors)
    if (length(bad_by) > 0L) {
      stop(
        "`contrast_by` names factor(s) not in the fit: ",
        paste(bad_by, collapse = ", "),
        call. = FALSE
      )
    }
  }

  .dd_validate_at(fit, at)
  grid <- .dd_build_emm_ref_grid(
    fit,
    at = at,
    factors_in_emm = factors_in_emm,
    validate = FALSE
  )
  emm_block <- .dd_brms_param_emms(
    fit,
    factors_in_emm = factors_in_emm,
    at = at,
    ci_level = ci_level
  )

  alpha2 <- (1 - ci_level) / 2
  ln10 <- log(10)
  empty_log10 <- tibble::tibble(
    contrast = character(),
    estimate = numeric(),
    std.error = numeric(),
    statistic = numeric(),
    df = numeric(),
    conf.low = numeric(),
    conf.high = numeric(),
    p.value = numeric(),
    post.prob = numeric()
  )
  empty_ratio <- tibble::tibble(
    contrast = character(),
    ratio = numeric(),
    conf.low = numeric(),
    conf.high = numeric(),
    p.value = numeric(),
    post.prob = numeric()
  )

  finish <- function(contrasts_log10, contrasts_ratio, by_applied) {
    block <- list(emmeans = emm_block, contrasts_log10 = contrasts_log10)
    if (report_ratios) {
      block$contrasts_ratio <- contrasts_ratio
    }
    results_list <- list(k = block)
    class(results_list) <- "beezdiscounting_comparison"
    attr(results_list, "backend") <- "brms"
    attr(results_list, "adjustment_method") <-
      "none (posterior summaries; see post.prob)"
    attr(results_list, "compare_specs_used") <- if (is.null(compare_specs)) {
      "all fitted factors"
    } else {
      deparse(compare_specs)
    }
    attr(results_list, "contrast_type_used") <- contrast_type
    attr(results_list, "contrast_by_used") <- if (
      is.null(contrast_by) ||
        !by_applied
    ) {
      "NULL"
    } else {
      paste(contrast_by, collapse = ", ")
    }
    results_list
  }

  if (isTRUE(grid$is_intercept_only) || nrow(grid$ref_X) < 2L) {
    return(finish(empty_log10, empty_ratio, FALSE))
  }

  draws <- .dd_brms_draws_matrix(fit)
  b <- as.matrix(
    draws[, .dd_brms_logk_draw_vars(fit, colnames(draws)), drop = FALSE]
  )
  ref_X <- grid$ref_X
  n <- nrow(ref_X)
  cell_lin <- b %*% t(ref_X)
  level_combos <- grid$level_combos
  use_factors <- grid$use_factors

  # contrast_by resolution mirroring the TMB path (Codex 041-B1): resolve
  # against the fitted factor set, abort loudly when a resolved by-var was
  # excluded by compare_specs, and ignore a redundant sole-factor by.
  effective_by <- character(0)
  if (!is.null(contrast_by)) {
    effective_by <- intersect(contrast_by, fitted_factors)
    if (length(effective_by) > 0L && !all(effective_by %in% use_factors)) {
      not_in <- setdiff(effective_by, use_factors)
      stop(
        "`contrast_by` factor(s) ",
        paste(not_in, collapse = ", "),
        " are not in `compare_specs`. Name the by-variable(s) in ",
        "`compare_specs` to condition contrasts on them.",
        call. = FALSE
      )
    }
    if (
      length(effective_by) > 0L &&
        length(use_factors) == 1L &&
        identical(sort(use_factors), sort(effective_by))
    ) {
      message(
        "  `contrast_by` (",
        paste(effective_by, collapse = ", "),
        ") is redundant with `compare_specs` for simple contrasts. ",
        "Ignoring `contrast_by`."
      )
      effective_by <- character(0)
    }
  }
  comparison_factors <- setdiff(use_factors, effective_by)
  label_f <- function(i, fs) .dd_brms_grid_label(grid, i, fs)

  if (length(effective_by) > 0L) {
    by_key <- do.call(
      paste,
      c(
        lapply(effective_by, function(f) as.character(level_combos[[f]])),
        list(sep = "\r")
      )
    )
    blocks <- lapply(unique(by_key), function(k) which(by_key == k))
  } else {
    blocks <- list(seq_len(n))
  }

  log10_parts <- list()
  ratio_parts <- list()
  for (rows in blocks) {
    m <- length(rows)
    if (m < 2L) {
      next
    }
    if (contrast_type == "pairwise") {
      cmb <- utils::combn(m, 2L)
      lhs <- rows[cmb[1L, ]]
      rhs <- rows[cmb[2L, ]]
    } else {
      lhs <- rows[seq.int(2L, m)]
      rhs <- rep(rows[1L], m - 1L)
    }
    for (k in seq_along(lhs)) {
      d <- as.numeric(cell_lin[, lhs[k]] - cell_lin[, rhs[k]])
      d10 <- d / ln10
      pd <- .dd_brms_post_prob(d)
      lab <- paste(
        label_f(lhs[k], comparison_factors),
        "-",
        label_f(rhs[k], comparison_factors)
      )
      l_row <- tibble::tibble(
        contrast = lab,
        estimate = stats::median(d10),
        std.error = stats::sd(d10),
        statistic = NA_real_,
        df = NA_real_,
        conf.low = unname(stats::quantile(d10, alpha2)),
        conf.high = unname(stats::quantile(d10, 1 - alpha2)),
        p.value = NA_real_,
        post.prob = pd
      )
      r_row <- tibble::tibble(
        contrast = lab,
        ratio = stats::median(exp(d)),
        conf.low = unname(stats::quantile(exp(d), alpha2)),
        conf.high = unname(stats::quantile(exp(d), 1 - alpha2)),
        p.value = NA_real_,
        post.prob = pd
      )
      if (length(effective_by) > 0L) {
        bc <- stats::setNames(
          lapply(effective_by, function(f) {
            as.character(level_combos[[f]][rows[1L]])
          }),
          effective_by
        )
        l_row <- dplyr::bind_cols(tibble::as_tibble(bc), l_row)
        r_row <- dplyr::bind_cols(tibble::as_tibble(bc), r_row)
      }
      log10_parts[[length(log10_parts) + 1L]] <- l_row
      ratio_parts[[length(ratio_parts) + 1L]] <- r_row
    }
  }

  if (length(log10_parts) == 0L) {
    return(finish(empty_log10, empty_ratio, FALSE))
  }
  finish(
    dplyr::bind_rows(log10_parts),
    dplyr::bind_rows(ratio_parts),
    length(effective_by) > 0L
  )
}

# --- print / summary ----------------------------------------------------------------

#' @export
print.beezdiscounting_brms <- function(x, digits = 4, ...) {
  cat("Bayesian Mixed-Effects Discounting Model (brms)\n")
  cat(strrep("=", 50), "\n")
  cat("Equation:", x$param_info$equation, " Family:", x$param_info$family, "\n")
  cat(
    "Chains:",
    x$mcmc_info$chains,
    " Iterations:",
    x$mcmc_info$iter,
    " (warmup",
    paste0(x$mcmc_info$warmup, ")"),
    "\n"
  )
  cat(
    "Subjects:",
    x$param_info$n_subjects,
    " Observations:",
    x$param_info$n_obs,
    "\n"
  )
  cat("Converged:", ifelse(isTRUE(x$converged), "Yes", "No"), "\n")
  cat("\nCoefficients (posterior medians, natural scale):\n")
  ct <- .dd_brms_coef_table(x, "natural")
  print(
    as.data.frame(ct[, c("term_display", "estimate", "std.error")]),
    digits = digits,
    row.names = FALSE
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
    digits = digits,
    row.names = FALSE
  )
  cat("\nVariance components:\n")
  print(x$variance_components, digits = digits, row.names = FALSE)
  invisible(x)
}
