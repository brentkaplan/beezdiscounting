# Estimated marginal means (and, later, comparisons) of the discount rate k.
#
# The discounting model is linear in the fixed-effect coefficients on the log
# scale: log k = X beta_k. Every marginal mean is therefore a deterministic
# averaging matrix A applied to a REBUILT design basis X_full, then exp()
# back-transformed -- no emmeans / emm_basis / recover_data machinery is needed.
# `%||%` is the single package-level definition in R/utils.R; do NOT redefine it
# here.

# Resolve user-requested retained factors against the fit's fitted factor set.
# Single-parameter analogue of beezdemand .tmb_resolve_retained_factors()
# (no collapse_levels / per-param suffixing in the discounting model).
.dd_resolve_retained_factors <- function(requested, fitted_factors) {
  unresolved <- setdiff(requested, fitted_factors)
  if (length(unresolved) > 0L) {
    cli::cli_abort(c(
      "{cli::qty(unresolved)}Requested factor{?s} {.val {unresolved}} {?is/are} not in the model.",
      "i" = "Fitted factors: {.val {fitted_factors}}."
    ))
  }
  unique(intersect(requested, fitted_factors))
}

# Validate the `at` list for the discounting EMM/comparison helpers. Aborts on
# unnamed entries, names not in (factors u continuous_covariates), factor values
# not observed, or non-finite continuous values; warns once on multi-value
# continuous entries (first value used). Single-parameter simplification of
# beezdemand .tmb_validate_at() (no param_scope / collapse alias logic).
.dd_validate_at <- function(fit_obj, at) {
  if (is.null(at)) return(invisible(NULL))
  if (is.null(names(at)) || any(!nzchar(names(at)))) {
    cli::cli_abort(
      "All elements of {.arg at} must be named (use {.code list(factor = level, cov = value)})."
    )
  }
  cov_names <- fit_obj$param_info$continuous_covariates %||% character(0)
  all_factors <- fit_obj$param_info$factors %||% character(0)
  all_factors <- all_factors[nzchar(all_factors) & !is.na(all_factors)]
  valid_names <- c(all_factors, cov_names)
  bad_names <- setdiff(names(at), valid_names)
  if (length(bad_names) > 0L) {
    cli::cli_abort(c(
      "Unknown name{?s} in {.arg at}: {.field {bad_names}}.",
      "i" = "Valid names are the fit's factors and continuous covariates: {.field {valid_names}}.",
      "x" = "Did you mistype a factor or covariate name?"
    ))
  }
  data_used <- fit_obj$data
  for (nm in names(at)) {
    v <- at[[nm]]
    if (length(v) < 1L) {
      cli::cli_abort(c(
        "{.field {nm}} has length 0.",
        "i" = "Each {.arg at} entry must be a non-empty vector."
      ))
    }
    if (nm %in% all_factors) {
      observed <- sort(unique(as.character(data_used[[nm]])))
      bad_vals <- setdiff(as.character(v), observed)
      if (length(bad_vals) > 0L) {
        cli::cli_abort(c(
          "{.field {nm}} = {.val {bad_vals}} not an observed level.",
          "i" = "Observed levels: {.val {observed}}."
        ))
      }
    } else {
      v_num <- suppressWarnings(as.numeric(v))
      if (any(is.na(v_num)) || any(!is.finite(v_num))) {
        cli::cli_abort(c(
          "{.field {nm}} value{?s} {.val {as.character(v)}} not finite numeric.",
          "i" = "Continuous-covariate {.arg at} entries must be a single finite numeric value."
        ))
      }
      if (length(v) > 1L) {
        cli::cli_warn(c(
          "{.arg at${nm}} has length {length(v)}; using first value {.val {v_num[1]}}.",
          "i" = "Pass a single numeric value per continuous covariate."
        ))
      }
    }
  }
  invisible(NULL)
}

# Build the conditioned reference grid (level_combos) and the averaging-matrix
# design (ref_X = A %*% X_full) for k EMMs/comparisons. Single-parameter port of
# beezdemand .tmb_build_emm_ref_grid() (R/tmb-methods.R:3186); the only fitted
# design is fit_obj$formula_details$X and the only factor set is
# fit_obj$param_info$factors.
.dd_build_emm_ref_grid <- function(
  fit_obj,
  at = NULL,
  factors_in_emm = NULL,
  validate = TRUE
) {
  cov_names <- fit_obj$param_info$continuous_covariates %||% character(0)
  fitted_factors <- fit_obj$param_info$factors %||% character(0)
  fitted_factors <- fitted_factors[nzchar(fitted_factors) & !is.na(fitted_factors)]

  use_factors <- fitted_factors
  if (!is.null(factors_in_emm)) {
    if (length(factors_in_emm) == 0L) {
      use_factors <- character(0)            # ~ 1: marginalize everything
    } else {
      use_factors <- .dd_resolve_retained_factors(factors_in_emm, fitted_factors)
    }
  }

  if (isTRUE(validate)) .dd_validate_at(fit_obj, at)

  retained_factors <- use_factors
  is_intercept_only <- length(fitted_factors) == 0L && length(cov_names) == 0L

  if (is_intercept_only) {
    return(list(
      level_combos = NULL, ref_X = NULL, use_factors = character(0),
      cov_names = character(0), is_intercept_only = TRUE
    ))
  }

  data_used <- fit_obj$data

  factor_level_set <- function(f) {
    lv <- levels(data_used[[f]])
    if (is.null(lv)) lv <- sort(unique(as.character(data_used[[f]])))
    if (!is.null(at) && f %in% names(at)) lv <- lv[lv %in% as.character(at[[f]])]
    lv
  }
  fitted_levels <- stats::setNames(
    lapply(fitted_factors, factor_level_set), fitted_factors
  )
  if (length(fitted_factors) > 0L &&
      any(vapply(fitted_levels, length, integer(1)) == 0L)) {
    cli::cli_abort(c(
      "{.arg at} filter produced an empty reference grid.",
      "i" = "Check that the supplied factor levels exist in the data."
    ))
  }

  make_key <- function(df, cols) {
    if (length(cols) == 0L) return(rep("", nrow(df)))
    do.call(paste, c(lapply(cols, function(cc) as.character(df[[cc]])),
                     list(sep = "\r")))
  }
  as_training_factor <- function(values, f) {
    factor(values, levels = levels(data_used[[f]]) %||%
             sort(unique(as.character(data_used[[f]]))))
  }

  # Full factorial grid over ALL fitted factors (equal-weight averaging target).
  if (length(fitted_factors) > 0L) {
    full_combos <- do.call(expand.grid, c(
      lapply(fitted_factors, function(f) as_training_factor(fitted_levels[[f]], f)),
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
    names(full_combos) <- fitted_factors
  } else {
    full_combos <- data_used[1L, integer(0), drop = FALSE]
  }

  # Retained reference grid: crossing of retained factors, ordered by level
  # index, filtered to OBSERVED combinations (semi_join analog).
  if (length(retained_factors) > 0L) {
    level_combos <- do.call(expand.grid, c(
      lapply(retained_factors, function(f) as_training_factor(fitted_levels[[f]], f)),
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
    names(level_combos) <- retained_factors
    ord <- do.call(order, lapply(retained_factors,
                                 function(f) as.integer(level_combos[[f]])))
    level_combos <- level_combos[ord, , drop = FALSE]
    observed_keys <- make_key(
      unique(data_used[, retained_factors, drop = FALSE]), retained_factors
    )
    level_combos <- level_combos[
      make_key(level_combos, retained_factors) %in% observed_keys, , drop = FALSE
    ]
    if (nrow(level_combos) == 0L) {
      cli::cli_abort(c(
        "{.arg at} filter produced an empty reference grid.",
        "i" = "Check that the supplied factor levels are not mutually exclusive."
      ))
    }
    rownames(level_combos) <- NULL
  } else {
    level_combos <- data_used[1L, integer(0), drop = FALSE]
  }

  # Continuous covariates: training mean unless overridden via `at`.
  if (length(cov_names) > 0L) {
    for (cv in cov_names) {
      cv_value <- mean(data_used[[cv]], na.rm = TRUE)
      if (!is.null(at) && cv %in% names(at)) cv_value <- as.numeric(at[[cv]][1])
      full_combos[[cv]] <- cv_value
      level_combos[[cv]] <- cv_value
    }
  }

  # Rebuild the grid basis through the SAME build_fixed_rhs route as the fit
  # (B1) and pin it to the STORED contrasts, then verify (and reorder to) the
  # fitted column set -- abort loudly on mismatch. Using the identical RHS +
  # contrasts route is what guarantees the EMM columns align with beta_k.
  fitted_X <- fit_obj$formula_details$X
  emm_rhs <- build_fixed_rhs(
    factors = fitted_factors,
    factor_interaction = fit_obj$param_info$factor_interaction,
    continuous_covariates = cov_names,
    data = data_used
  )
  X_full <- stats::model.matrix(
    emm_rhs,
    data = full_combos,
    contrasts.arg = fit_obj$formula_details$contrasts
  )
  fitted_cols <- colnames(fitted_X)
  if (!is.null(fitted_cols)) {
    if (!setequal(colnames(X_full), fitted_cols)) {
      cli::cli_abort(c(
        "Could not reproduce the fitted design matrix for the EMM grid.",
        "i" = "Rebuilt columns: {.val {colnames(X_full)}}.",
        "i" = "Fitted columns: {.val {fitted_cols}}.",
        "x" = "This can happen if the model's factor levels or contrasts changed after fitting."
      ))
    }
    X_full <- X_full[, fitted_cols, drop = FALSE]
  }

  # Averaging matrix A (n_retained x n_full): equal weight 1/m on the m
  # full-grid rows matching each retained cell. ref_X = A %*% X_full keeps
  # ncol == length(beta_k).
  full_keys <- make_key(full_combos, retained_factors)
  ret_keys <- make_key(level_combos, retained_factors)
  A <- matrix(0, nrow = nrow(level_combos), ncol = nrow(full_combos))
  for (r in seq_len(nrow(level_combos))) {
    sel <- which(full_keys == ret_keys[r])
    A[r, sel] <- 1 / length(sel)
  }
  ref_X <- A %*% X_full
  colnames(ref_X) <- colnames(X_full)

  list(
    level_combos = level_combos,
    ref_X = ref_X,
    use_factors = retained_factors,
    cov_names = cov_names,
    is_intercept_only = FALSE
  )
}


# Resolve the beta_k covariance block, honoring the se_available gate (B2).
# - se_available FALSE (non-PD Hessian or sdreport failed): return an NA-filled
#   matrix and warn once, so EMM/contrast SE/CI/statistic/p come out NA rather
#   than unreliable Wald values.
# - se_available TRUE with sdr$cov.fixed present: return the beta_k block.
# - se_available TRUE but cov.fixed missing (rare): diagonal fallback, with a
#   warning that off-diagonal covariance is lost -- the EMM SEs stay valid but
#   contrast SEs (dx^T V dx) would be wrong.
.dd_resolve_beta_vcov <- function(fit, beta) {
  p <- length(beta)
  if (!isTRUE(fit$se_available)) {
    cli::cli_warn(c(
      "!" = "Standard errors are unreliable (non-PD Hessian or sdreport \\
             failed); reporting estimates with {.val NA} uncertainty.",
      "i" = "Interpret the point estimates only; check convergence / identifiability."
    ))
    return(matrix(NA_real_, p, p))
  }
  sdr <- fit$sdr
  if (!is.null(sdr) && !is.null(sdr$cov.fixed)) {
    full_vcov <- as.matrix(sdr$cov.fixed)
    target_idx <- which(names(fit$opt$par) == "beta_k")
    if (length(target_idx) == p) {
      return(full_vcov[target_idx, target_idx, drop = FALSE])
    }
  }
  cli::cli_warn(
    "Using a diagonal covariance fallback; EMM and contrast uncertainty may be \\
     approximate / unreliable (off-diagonal covariance unavailable)."
  )
  se_vals <- unname(fit$model$se[names(fit$model$coefficients) == "beta_k"])
  diag(se_vals^2, nrow = p)
}


#' Estimated marginal means of the discount rate \code{k}
#'
#' @description
#' Computes estimated marginal means (EMMs) of the discount rate \code{k} from a
#' fitted \code{beezdiscounting_tmb} model. EMMs are computed on the
#' \code{log k} scale (linear in the fixed-effect coefficients) using the
#' averaging-matrix reference grid, then back-transformed with \code{exp()} so
#' that \code{k = exp(k_log)}. Standard errors use the \code{beta_k} block of
#' \code{TMB::sdreport()}'s fixed-effect covariance; intervals are Wald on the
#' log scale and exponentiated.
#'
#' @param fit A \code{beezdiscounting_tmb} object.
#' @param factors_in_emm Character vector of factors to retain in the EMM
#'   reference grid. A strict subset marginalizes the omitted factors with equal
#'   weights across the full crossing of their levels (emmeans' default
#'   \code{weights = "equal"}); \code{NULL} (default) retains all fitted factors;
#'   \code{character(0)} marginalizes everything to a single grand-mean cell.
#' @param at Named list specifying factor levels and/or continuous-covariate
#'   values for conditional EMMs (one numeric per covariate; multiple values
#'   warn and use the first). \code{at} on an omitted factor restricts the
#'   level set averaged over.
#' @param ci_level Numeric confidence level for intervals (default 0.95).
#' @param ... Additional arguments (currently unused).
#'
#' @return A tibble with columns \code{level}, \code{k}, \code{k_log},
#'   \code{std.error}, \code{conf.low}, \code{conf.high}. \code{k_log} is the
#'   marginal mean on the \code{log k} scale; \code{k = exp(k_log)};
#'   \code{std.error} is the SE of \code{k_log}; the intervals are on the
#'   \code{k} (natural) scale.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dd_ip(
#'   n_subjects = 30, n_conditions = 2, delta_k = c(0, log(3)), seed = 1
#' )
#' fit <- fit_dd_tmb(sim, factors = "condition", verbose = 0)
#' get_dd_param_emms(fit)
#' }
#'
#' @export
get_dd_param_emms <- function(
  fit,
  factors_in_emm = NULL,
  at = NULL,
  ci_level = 0.95,
  ...
) {
  if (inherits(fit, "beezdiscounting_choice") &&
      identical(fit$param_info$mode, "descriptive")) {
    cli::cli_abort(c(
      "emmeans is not available for descriptive (Young 2018) choice models.",
      "i" = "The descriptive model has no discount rate {.code k}; use \\
             {.fn VarCorr} / {.fn ranef} for per-subject slope (co)variances."))
  }
  coefs <- fit$model$coefficients
  beta_idx <- which(names(coefs) == "beta_k")
  beta <- unname(coefs[beta_idx])

  # beta_k covariance block, gated on se_available (B2): NA uncertainty when
  # SEs are unreliable (non-PD Hessian / sdreport failed), the sdreport block
  # otherwise.
  vcov_mat <- .dd_resolve_beta_vcov(fit, beta)

  .dd_validate_at(fit, at)

  grid <- .dd_build_emm_ref_grid(
    fit, at = at, factors_in_emm = factors_in_emm, validate = FALSE
  )
  z <- stats::qnorm((1 + ci_level) / 2)

  if (isTRUE(grid$is_intercept_only)) {
    est <- beta[1L]
    se <- sqrt(vcov_mat[1L, 1L])
    return(tibble::tibble(
      level = "(Intercept)",
      k = exp(est),
      k_log = est,
      std.error = se,
      conf.low = exp(est - z * se),
      conf.high = exp(est + z * se)
    ))
  }

  use_factors <- grid$use_factors
  cov_names <- grid$cov_names
  level_combos <- grid$level_combos
  ref_X <- grid$ref_X

  if (ncol(ref_X) != length(beta)) {
    cli::cli_abort(c(
      "Reference-grid design has {ncol(ref_X)} column{?s} but the fitted \\
       coefficient vector has {length(beta)}.",
      "x" = "Design basis mismatch; cannot evaluate EMMs."
    ))
  }

  level_label_for <- function(i) {
    if (length(use_factors) > 0L) {
      paste(vapply(use_factors, function(f)
        paste0(f, "=", level_combos[[f]][i]), character(1)), collapse = ", ")
    } else if (length(cov_names) > 0L) {
      paste(vapply(cov_names, function(cv)
        paste0(cv, "=", level_combos[[cv]][i]), character(1)), collapse = ", ")
    } else {
      "(Intercept)"
    }
  }

  cell_est <- as.numeric(ref_X %*% beta)
  cell_se <- sqrt(diag(ref_X %*% vcov_mat %*% t(ref_X)))

  tibble::tibble(
    level = vapply(seq_len(nrow(ref_X)), level_label_for, character(1)),
    k = exp(cell_est),
    k_log = cell_est,
    std.error = cell_se,
    conf.low = exp(cell_est - z * cell_se),
    conf.high = exp(cell_est + z * cell_se)
  )
}


#' Factor-level comparisons of the discount rate \code{k}
#'
#' @description
#' Computes factor-level contrasts of the discount rate \code{k} from a fitted
#' \code{beezdiscounting_tmb} model. \code{k} is linear in the fixed-effect
#' coefficients on the natural-log scale (\code{log k = X beta_k}), so each
#' contrast is a linear combination of \code{beta_k} with a Wald standard error
#' from the \code{beta_k} block of \code{TMB::sdreport()}'s fixed-effect
#' covariance. Contrasts are \strong{reported on the \code{log10} scale} and,
#' optionally, as multiplicative ratios (\code{ratio = exp(est_log)}). The
#' returned container mirrors the beezdemand \code{beezdemand_comparison} shape,
#' so [tidy.beezdiscounting_comparison()] gives a flat, cross-backend frame.
#'
#' @param fit A \code{beezdiscounting_tmb} object.
#' @param compare_specs Optional one-sided formula naming the factor subset to
#'   contrast (e.g. \code{~ condition}). Omitted fitted factors are marginalized
#'   over with equal weights across the full crossing of their levels. If
#'   \code{NULL} (default), all fitted factors are retained. Unknown names abort.
#' @param contrast_type Character. \code{"pairwise"} (all \code{choose(n, 2)}
#'   pairs, factor-level order) or \code{"trt.vs.ctrl"} (each level vs. the
#'   first/reference level).
#' @param contrast_by Optional \code{NULL} (default) or character vector of
#'   factor name(s) within \code{compare_specs} to condition the contrasts on.
#'   Within each observed combination of by-level(s), pairwise (or
#'   \code{trt.vs.ctrl}) contrasts are computed over the remaining (non-by)
#'   factors, with p-value adjustment applied \strong{per by-cell}. A
#'   \code{contrast_by} factor absent from \code{compare_specs} aborts.
#' @param adjust Character. P-value adjustment method; must be one of
#'   \code{stats::p.adjust.methods} (default \code{"holm"}). emmeans-only methods
#'   such as \code{"tukey"}/\code{"sidak"} are rejected (the TMB backend uses an
#'   asymptotic \emph{z} + \code{stats::p.adjust()}).
#' @param at Named list specifying factor levels and/or continuous-covariate
#'   values to condition on, as in [get_dd_param_emms()].
#' @param ci_level Numeric. Confidence level for intervals (default 0.95).
#' @param report_ratios Logical. If \code{TRUE} (default), include a
#'   \code{contrasts_ratio} block (multiplicative ratios).
#' @param ... Additional arguments (reserved; \code{factors_in_emm} is accepted
#'   as a lower-level alternative to \code{compare_specs}).
#'
#' @return A \code{beezdiscounting_comparison} object: a list with a single
#'   element \code{k}, itself a list with \code{emmeans} (the EMM tibble from
#'   [get_dd_param_emms()]), \code{contrasts_log10} (columns \code{contrast},
#'   \code{estimate}, \code{std.error}, \code{statistic}, \code{df},
#'   \code{conf.low}, \code{conf.high}, \code{p.value}), and (if
#'   \code{report_ratios}) \code{contrasts_ratio} (columns \code{contrast},
#'   \code{ratio}, \code{conf.low}, \code{conf.high}, \code{p.value}). When
#'   \code{contrast_by} is active, the contrast tables gain leading by-column(s)
#'   before \code{contrast}. Attributes \code{backend}, \code{adjustment_method},
#'   \code{compare_specs_used}, \code{contrast_type_used}, and
#'   \code{contrast_by_used} describe the call.
#'
#' @seealso [tidy.beezdiscounting_comparison()] for the cross-backend frame.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dd_ip(
#'   n_subjects = 30, n_conditions = 2, delta_k = c(0, log(3)), seed = 1
#' )
#' fit <- fit_dd_tmb(sim, factors = "condition", verbose = 0)
#' get_dd_comparisons(fit, contrast_type = "pairwise")
#' }
#'
#' @export
get_dd_comparisons <- function(
  fit,
  compare_specs = NULL,
  contrast_type = c("pairwise", "trt.vs.ctrl"),
  contrast_by = NULL,
  adjust = "holm",
  at = NULL,
  ci_level = 0.95,
  report_ratios = TRUE,
  ...
) {
  if (inherits(fit, "beezdiscounting_choice") &&
      identical(fit$param_info$mode, "descriptive")) {
    cli::cli_abort(c(
      "emmeans is not available for descriptive (Young 2018) choice models.",
      "i" = "The descriptive model has no discount rate {.code k}; use \\
             {.fn VarCorr} / {.fn ranef} for per-subject slope (co)variances."))
  }
  contrast_type <- match.arg(contrast_type)

  fitted_factors <- fit$param_info$factors %||% character(0)
  fitted_factors <- fitted_factors[nzchar(fitted_factors) & !is.na(fitted_factors)]

  # adjust: validate against the base-R set. emmeans-only methods
  # (tukey/sidak/scheffe/mvt) are not implementable with stats::p.adjust().
  if (!isTRUE(adjust %in% stats::p.adjust.methods)) {
    cli::cli_abort(c(
      "{.arg adjust} = {.val {adjust}} is not a valid p-value adjustment method.",
      "i" = "Valid methods: {.val {stats::p.adjust.methods}}.",
      "x" = "emmeans-only methods (e.g. {.val tukey}, {.val sidak}) are \\
             unavailable on the TMB backend (asymptotic z + \\
             {.fn stats::p.adjust})."
    ))
  }

  # Resolve the retained factor set: `compare_specs` (canonical) wins, else the
  # lower-level `factors_in_emm` via `...` (backward compatible).
  dots <- list(...)
  factors_in_emm <- NULL
  if (!is.null(compare_specs)) {
    if (!inherits(compare_specs, "formula")) {
      cli::cli_abort("{.arg compare_specs} must be a one-sided formula (e.g. {.code ~ condition}).")
    }
    factors_in_emm <- all.vars(compare_specs)
    bad <- setdiff(factors_in_emm, fitted_factors)
    if (length(bad) > 0L) {
      cli::cli_abort(c(
        "{.arg compare_specs} names factor{?s} not in the fit: {.val {bad}}.",
        "i" = "Fitted factors: {.val {fitted_factors}}."
      ))
    }
  } else if (!is.null(dots$factors_in_emm)) {
    factors_in_emm <- dots$factors_in_emm
  }

  # contrast_by: NULL or character(1+). Boundary validation (loud) catches typos
  # against the fitted factor set; resolution against `compare_specs` happens in
  # the worker (single-parameter model: no collapse aliasing).
  if (!is.null(contrast_by)) {
    if (!is.character(contrast_by)) {
      cli::cli_abort(
        "{.arg contrast_by} must be {.code NULL} or a character vector of \\
         factor name(s)."
      )
    }
    if (length(contrast_by) == 0L) {
      contrast_by <- NULL
    } else {
      bad_by <- setdiff(contrast_by, fitted_factors)
      if (length(bad_by) > 0L) {
        cli::cli_abort(c(
          "{.arg contrast_by} names factor{?s} not in the fit: {.val {bad_by}}.",
          "i" = "Fitted factors: {.val {fitted_factors}}."
        ))
      }
    }
  }

  # Validate `at` once at the public boundary (single multi-value warning).
  .dd_validate_at(fit, at)

  block <- .dd_compare_k(
    fit, factors_in_emm, contrast_type, adjust, at, ci_level,
    report_ratios, contrast_by
  )

  results_list <- list(k = block)
  class(results_list) <- "beezdiscounting_comparison"
  attr(results_list, "backend") <- "tmb"
  attr(results_list, "adjustment_method") <- adjust
  attr(results_list, "compare_specs_used") <- if (is.null(compare_specs)) {
    "all fitted factors"
  } else {
    deparse(compare_specs)
  }
  attr(results_list, "contrast_type_used") <- contrast_type
  # Record the user-requested name(s) only when by-grouping was actually applied
  # (otherwise "NULL" so the flattener/print do not synthesize an all-NA col).
  by_applied <- isTRUE(attr(block, "by_applied"))
  attr(results_list, "contrast_by_used") <- if (is.null(contrast_by) || !by_applied) {
    "NULL"
  } else {
    paste(contrast_by, collapse = ", ")
  }
  attr(results_list$k, "by_applied") <- NULL
  results_list
}

# Build the single-parameter (k) comparison block for the TMB discounting
# backend. Returns list(emmeans, contrasts_log10[, contrasts_ratio]); carries a
# `by_applied` attribute so the public function can decide whether to surface
# the `contrast_by_used` metadata. The discount rate is linear in beta_k on the
# natural-log scale; contrasts are reported on the log10 scale and (optionally)
# as multiplicative ratios on the natural scale.
.dd_compare_k <- function(fit, factors_in_emm, contrast_type, adjust, at,
                          ci_level, report_ratios, contrast_by = NULL) {
  coefs <- fit$model$coefficients
  beta <- unname(coefs[names(coefs) == "beta_k"])

  # beta_k covariance block, gated on se_available (B2).
  vcov_mat <- .dd_resolve_beta_vcov(fit, beta)

  grid <- .dd_build_emm_ref_grid(
    fit, at = at, factors_in_emm = factors_in_emm, validate = FALSE
  )
  z <- stats::qnorm((1 + ci_level) / 2)
  ln10 <- log(10)

  # The native EMM block (same tibble get_dd_param_emms returns). `at` was
  # already validated (and any multi-value warning emitted once) at the public
  # boundary, so suppress the inner re-validation's duplicate warning.
  emm_block <- withCallingHandlers(
    get_dd_param_emms(
      fit, factors_in_emm = factors_in_emm, at = at, ci_level = ci_level
    ),
    warning = function(w) {
      # Muffle the inner EMM's duplicate warnings: the multi-value `at` notice
      # (already emitted once at the public boundary) and the B2 se-unreliable
      # warning (this worker emits its own via .dd_resolve_beta_vcov above).
      if (grepl("using first value|unreliable", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )

  empty_log10 <- tibble::tibble(
    contrast = character(), estimate = numeric(), std.error = numeric(),
    statistic = numeric(), df = numeric(),
    conf.low = numeric(), conf.high = numeric(), p.value = numeric()
  )
  empty_ratio <- tibble::tibble(
    contrast = character(), ratio = numeric(),
    conf.low = numeric(), conf.high = numeric(), p.value = numeric()
  )
  finish_empty <- function() {
    out <- list(emmeans = emm_block, contrasts_log10 = empty_log10)
    if (report_ratios) out$contrasts_ratio <- empty_ratio
    attr(out, "by_applied") <- FALSE
    out
  }

  # An intercept-only grid (or any single-cell grid) admits no contrasts.
  if (isTRUE(grid$is_intercept_only) || nrow(emm_block) < 2L) {
    return(finish_empty())
  }

  use_factors <- grid$use_factors
  cov_names <- grid$cov_names
  level_combos <- grid$level_combos
  ref_X <- grid$ref_X
  n <- nrow(ref_X)

  # Cross-backend contrast labels are built from STRUCTURED ref-grid level
  # values, parametrized by the factor subset used (by-vars are excluded so the
  # within-cell label matches the `at = `-filtered route).
  native_label_f <- function(i, fs) {
    if (length(fs) > 0L) {
      paste(vapply(fs, function(f)
        paste0(f, "=", as.character(level_combos[[f]][i])), character(1)),
        collapse = ", ")
    } else if (length(cov_names) > 0L) {
      paste(vapply(cov_names, function(cv)
        paste0(cv, "=", level_combos[[cv]][i]), character(1)), collapse = ", ")
    } else {
      "(Intercept)"
    }
  }

  # Resolve contrast_by against this model's retained factor set. Single-
  # parameter discounting model: no collapse aliasing, so a by-var is either in
  # `use_factors` or it aborts (it was already validated against the fit's
  # factors at the public boundary).
  effective_by <- character(0)
  if (n >= 2L && !is.null(contrast_by)) {
    if (!all(contrast_by %in% use_factors)) {
      not_in <- setdiff(contrast_by, use_factors)
      cli::cli_abort(c(
        "{cli::qty(not_in)}{.arg contrast_by} factor{?s} {.val {not_in}} \\
         {?is/are} not in {.arg compare_specs}.",
        "i" = "{cli::qty(use_factors)}{.arg compare_specs} factor{?s}: {.val {use_factors}}.",
        "x" = "Name the by-variable(s) in {.arg compare_specs} to condition contrasts on them."
      ))
    }
    effective_by <- contrast_by
    # Redundant-by: a length-1 compare_specs equal to the by-set yields no
    # remaining comparison factor; ignore the by-grouping for this case.
    if (length(use_factors) == 1L &&
        identical(sort(use_factors), sort(effective_by))) {
      cli::cli_inform(c(
        "i" = "{.arg contrast_by} = {.val {contrast_by}} matches the sole \\
               comparison factor and was ignored.",
        " " = "The contrasts are computed over all levels of \\
               {.val {use_factors}} without by-grouping."
      ))
      effective_by <- character(0)
    }
  }

  comparison_factors <- setdiff(use_factors, effective_by)

  # Row blocks: one per observed by-cell (factor-level order preserved from the
  # ref grid), or a single global block when no by-grouping is active.
  if (length(effective_by) > 0L) {
    by_key <- do.call(paste, c(
      lapply(effective_by, function(f) as.character(level_combos[[f]])),
      list(sep = "\r")
    ))
    blocks <- lapply(unique(by_key), function(k) which(by_key == k))
  } else {
    blocks <- list(seq_len(n))
  }

  # Compute pairwise / trt.vs.ctrl contrasts within one block of grid rows.
  # p-values are adjusted WITHIN the block (per by-cell).
  do_block <- function(rows) {
    m <- length(rows)
    if (m < 2L) return(NULL)
    if (contrast_type == "pairwise") {
      cmb <- utils::combn(m, 2L)
      lhs <- rows[cmb[1L, ]]
      rhs <- rows[cmb[2L, ]]
    } else {
      lhs <- rows[seq.int(2L, m)]
      rhs <- rep(rows[1L], m - 1L)
    }
    est_log <- numeric(length(lhs))
    se_log <- numeric(length(lhs))
    label <- character(length(lhs))
    for (k in seq_along(lhs)) {
      dx <- ref_X[lhs[k], ] - ref_X[rhs[k], ]
      est_log[k] <- sum(dx * beta)
      se_log[k] <- sqrt(as.numeric(t(dx) %*% vcov_mat %*% dx))
      label[k] <- paste(native_label_f(lhs[k], comparison_factors), "-",
                        native_label_f(rhs[k], comparison_factors))
    }
    zstat <- est_log / se_log
    p_adj <- stats::p.adjust(2 * stats::pnorm(-abs(zstat)), method = adjust)
    est_log10 <- est_log / ln10
    se_log10 <- se_log / ln10
    list(
      log10 = tibble::tibble(
        contrast = label, estimate = est_log10, std.error = se_log10,
        statistic = zstat, df = Inf,
        conf.low = est_log10 - z * se_log10,
        conf.high = est_log10 + z * se_log10, p.value = p_adj
      ),
      ratio = tibble::tibble(
        contrast = label, ratio = exp(est_log),
        conf.low = exp(est_log - z * se_log),
        conf.high = exp(est_log + z * se_log), p.value = p_adj
      ),
      first_row = rows[1L]
    )
  }

  # Build the by-column tibble for a block (user-requested by names).
  by_cols_for <- function(first_row, nrows) {
    if (length(effective_by) == 0L) return(NULL)
    cols <- lapply(effective_by, function(f) {
      rep(as.character(level_combos[[f]][first_row]), nrows)
    })
    tibble::as_tibble(stats::setNames(cols, effective_by))
  }

  block_results <- Filter(Negate(is.null), lapply(blocks, do_block))

  if (length(block_results) == 0L) {
    res <- finish_empty()
    attr(res, "by_applied") <- length(effective_by) > 0L
    return(res)
  }

  log10_parts <- lapply(block_results, function(r) {
    bc <- by_cols_for(r$first_row, nrow(r$log10))
    if (is.null(bc)) r$log10 else dplyr::bind_cols(bc, r$log10)
  })
  ratio_parts <- lapply(block_results, function(r) {
    bc <- by_cols_for(r$first_row, nrow(r$ratio))
    if (is.null(bc)) r$ratio else dplyr::bind_cols(bc, r$ratio)
  })

  out <- list(
    emmeans = emm_block,
    contrasts_log10 = dplyr::bind_rows(log10_parts)
  )
  if (report_ratios) out$contrasts_ratio <- dplyr::bind_rows(ratio_parts)
  attr(out, "by_applied") <- length(effective_by) > 0L
  out
}


#' Tidy a discounting comparison into a flat contrasts frame
#'
#' @description
#' [generics::tidy()] method for \code{beezdiscounting_comparison} objects
#' (returned by [get_dd_comparisons()]). Produces a flat long tibble whose
#' \strong{column names and order} match the beezdemand
#' \code{beezdemand_comparison} tidier, enabling downstream consumers to bind
#' rows from both backends into a single frame. Note that the \code{contrast}
#' label dialect differs: beezdiscounting emits native \code{factor=level}
#' labels (e.g. \code{"condition=C1 - condition=C2"}), while beezdemand emits
#' bare level labels (e.g. \code{"C1 - C2"}). Fully uniform emmeans-style
#' contrast labels across backends are a future cross-backend item.
#' The nested object keeps the native dialect (see [get_dd_comparisons()]).
#'
#' @param x A \code{beezdiscounting_comparison} object.
#' @param exponentiate Logical. If \code{TRUE}, return base-invariant ratios
#'   (\code{estimate = 10^estimate}, CIs back-transformed); \code{std.error}
#'   becomes \code{NA} following broom's convention for exponentiated fits.
#'   Default \code{FALSE} (log10-scale contrasts).
#' @param ... Unused.
#'
#' @return A tibble with columns \code{param}, \code{contrast},
#'   \code{estimate}, \code{std.error}, \code{statistic}, \code{df},
#'   \code{conf.low}, \code{conf.high}, \code{p.value} (with leading by-column(s)
#'   inserted before \code{param} when \code{contrast_by} is active). Estimates
#'   and CIs are on the log10 scale (or ratios when \code{exponentiate = TRUE}).
#'   \code{statistic} is an asymptotic \emph{z} (\code{df = Inf}).
#'
#' @importFrom generics tidy
#' @export
tidy.beezdiscounting_comparison <- function(x, exponentiate = FALSE, ...) {
  by_used <- attr(x, "contrast_by_used")
  by_active <- !(is.null(by_used) || length(by_used) == 0L ||
                   identical(by_used, "") || identical(by_used, "NULL"))
  by_names <- if (by_active) trimws(strsplit(by_used, ",")[[1]]) else character(0)

  base_cols <- list(
    param = character(), contrast = character(), estimate = numeric(),
    std.error = numeric(), statistic = numeric(), df = numeric(),
    conf.low = numeric(), conf.high = numeric(), p.value = numeric()
  )

  rows <- lapply(names(x), function(p) {
    cl <- x[[p]]$contrasts_log10
    if (is.null(cl) || nrow(cl) == 0L || !("estimate" %in% names(cl))) {
      return(NULL)
    }
    base <- tibble::tibble(
      param = p, contrast = cl$contrast,
      estimate = cl$estimate, std.error = cl$std.error,
      statistic = cl$statistic, df = cl$df,
      conf.low = cl$conf.low, conf.high = cl$conf.high,
      p.value = cl$p.value
    )
    if (by_active) {
      by_cols <- lapply(by_names, function(nm) {
        if (nm %in% names(cl)) as.character(cl[[nm]]) else rep(NA_character_, nrow(cl))
      })
      base <- dplyr::bind_cols(
        tibble::as_tibble(stats::setNames(by_cols, by_names)), base
      )
    }
    base
  })

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0L) {
    out <- tibble::as_tibble(base_cols)
    if (by_active) {
      by_empty <- stats::setNames(
        rep(list(character()), length(by_names)), by_names
      )
      out <- dplyr::bind_cols(tibble::as_tibble(by_empty), out)
    }
  }

  if (isTRUE(exponentiate)) {
    # Base-invariant ratios; std.error is NA per broom's exponentiated-fit
    # convention (the delta-method SE does not transform multiplicatively).
    out$estimate <- 10^out$estimate
    out$conf.low <- 10^out$conf.low
    out$conf.high <- 10^out$conf.high
    out$std.error <- NA_real_
  }
  out
}
