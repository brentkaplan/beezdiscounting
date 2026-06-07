#' @useDynLib beezdiscounting, .registration = TRUE
NULL

# ==============================================================================
# Internal Helper Functions for TMB Mixed-Effects Discounting Models
# ==============================================================================

#' Prepare data for the TMB mixed-effects discounting model
#'
#' Operates on the FULL model frame from `.dd_validate_ip()` (canonical
#' `id`/`x`/`y` PLUS any retained factor/covariate columns). It complete-cases
#' that frame ONCE over all modeling columns (`id`, `x`, `y`, and `extra_cols`),
#' then derives the parallel arrays `y`/`x`/`subject_id` from the same filtered
#' rows so the later design matrix, TMB arrays, and `fit$data` all share one
#' row order (row-coherence; mirrors beezdemand's `data_for_design`). Builds a
#' 0-indexed `subject_id` aligned to `subject_levels` (the C++ template indexes
#' `u(subject_id, 0)` from 0).
#'
#' @param data Data frame already coerced/clamped by `.dd_validate_ip()`,
#'   retaining factor/covariate columns alongside id/x/y.
#' @param y_var,x_var,id_var Character column names.
#' @param extra_cols Character vector of factor/covariate column names that must
#'   also be free of NA for a row to be kept (the union of `factors` and
#'   `continuous_covariates`). Defaults to none.
#' @return A list with `y`, `x`, `subject_id` (0-indexed integer),
#'   `subject_levels`, `n_subjects`, `n_obs`, and `data` (the SINGLE complete-
#'   cased model frame: canonical `id`, `x`, `y` first, then the retained
#'   factor/covariate columns, in `data`'s row order).
#' @keywords internal
.dd_tmb_prepare_data <- function(data, y_var = "y", x_var = "x", id_var = "id",
                                 extra_cols = NULL) {
  # Canonicalize id/x/y in place, keeping every other column intact.
  frame <- data
  frame[[id_var]] <- as.character(frame[[id_var]])
  frame[[x_var]] <- as.numeric(frame[[x_var]])
  frame[[y_var]] <- as.numeric(frame[[y_var]])

  # Complete-case ONCE over all modeling columns (id/x/y + extras). Doing this
  # in a single place is load-bearing: model.matrix() must NOT drop a different
  # set of rows than the prepared arrays.
  model_cols <- unique(c(id_var, x_var, y_var, extra_cols))
  model_cols <- intersect(model_cols, names(frame))
  keep <- stats::complete.cases(frame[, model_cols, drop = FALSE])
  frame <- frame[keep, , drop = FALSE]

  if (nrow(frame) == 0L) {
    stop("No complete cases remain after dropping NA rows.", call. = FALSE)
  }

  ids <- frame[[id_var]]
  x <- frame[[x_var]]
  y <- frame[[y_var]]

  # 0-indexed subject map, levels sorted for reproducibility
  subject_levels <- sort(unique(ids))
  n_subjects <- length(subject_levels)
  subject_map <- stats::setNames(seq_along(subject_levels) - 1L, subject_levels)
  subject_id <- as.integer(subject_map[ids])

  # Single filtered frame: canonical id/x/y first, then retained extras, in the
  # surviving row order. Drop unused factor levels so downstream model.matrix /
  # ref grids see only present levels.
  extras <- setdiff(intersect(extra_cols, names(frame)), c(id_var, x_var, y_var))
  cleaned <- data.frame(
    id = ids,
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
  for (col in extras) {
    v <- frame[[col]]
    if (is.factor(v)) v <- droplevels(v)
    cleaned[[col]] <- v
  }

  list(
    y = y,
    x = x,
    subject_id = subject_id,
    subject_levels = subject_levels,
    n_subjects = n_subjects,
    n_obs = length(y),
    data = cleaned
  )
}


#' Build the fixed-effect design matrix for log k
#'
#' Constructs the `model.matrix` for the `log k = Xbeta` linear predictor from
#' between-subject factors, an optional pairwise interaction, and continuous
#' covariates. With no inputs the design is intercept-only.
#'
#' @param data Cleaned long data frame (canonical `id`/`x`/`y` plus any factor
#'   or covariate columns).
#' @param factors Character vector of factor column names, or `NULL`.
#' @param factor_interaction Logical; if `TRUE` and >= 2 factors, include their
#'   interaction (uses `*`); otherwise main effects (`+`).
#' @param continuous_covariates Character vector of covariate names, or `NULL`.
#' @return A list with `X` (the model matrix for log-k FE), `rhs` (the one-sided
#'   RHS **formula** from [build_fixed_rhs()]), and `contrasts`
#'   (`attr(X, "contrasts")`, the per-factor contrast scheme). The `rhs` and
#'   `contrasts` are stored so `.dd_build_emm_ref_grid()` can rebuild a
#'   column-aligned design via the same route.
#' @keywords internal
.dd_tmb_build_design <- function(data, factors = NULL,
                                 factor_interaction = FALSE,
                                 continuous_covariates = NULL) {
  # Single source of the RHS (shared with the emmeans reference grid in E.1).
  rhs <- build_fixed_rhs(
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    data = data
  )
  X <- stats::model.matrix(rhs, data = data)

  list(X = X, rhs = rhs, contrasts = attr(X, "contrasts"))
}


#' Build the TMB data list for MixedDiscounting
#'
#' @param prepared Output from [.dd_tmb_prepare_data()].
#' @param design Output from [.dd_tmb_build_design()].
#' @param equation One of "mazur", "exponential".
#' @param family One of "sltb", "gaussian".
#' @return A list whose names match the C++ `DATA_*` macros: `model`, `y`, `x`,
#'   `subject_id` (0-indexed integer), `X`, `eqn_type`, `family`, `n_obs`,
#'   `n_subjects`.
#' @keywords internal
.dd_tmb_build_tmb_data <- function(prepared, design, equation, family) {
  eqn_type <- switch(equation,
    mazur = 0L,
    exponential = 1L,
    stop("Unknown equation: ", equation, call. = FALSE)
  )
  fam_type <- switch(family,
    sltb = 0L,
    gaussian = 1L,
    stop("Unknown family: ", family, call. = FALSE)
  )

  list(
    model = "MixedDiscounting",
    y = as.numeric(prepared$y),
    x = as.numeric(prepared$x),
    subject_id = as.integer(prepared$subject_id),
    X = as.matrix(design$X),
    eqn_type = eqn_type,
    family = fam_type,
    n_obs = as.integer(prepared$n_obs),
    n_subjects = as.integer(prepared$n_subjects)
  )
}


#' Generate default starting values for the TMB discounting model
#'
#' `beta_k` intercept is data-driven: invert the discounting function at the
#' median `y` observed at the minimum delay. `log_sigma_u = log(0.5)`. The
#' generic auxiliary scalar `log_aux` starts at `log(8)` for SLT-beta
#' (population precision phi) and `log(0.1)` for Gaussian (residual SD).
#'
#' @param prepared Output from [.dd_tmb_prepare_data()].
#' @param design Output from [.dd_tmb_build_design()].
#' @param family One of "sltb", "gaussian".
#' @param equation One of "mazur", "exponential" (for the start inversion).
#' @return A parameters list: `beta_k` (length `ncol(X)`), `log_sigma_u`,
#'   `log_aux`, `u` (matrix `n_subjects` x 1).
#' @keywords internal
.dd_tmb_default_starts <- function(prepared, design, family,
                                   equation = "mazur") {
  y <- prepared$y
  x <- prepared$x
  p <- ncol(design$X)
  n_subjects <- prepared$n_subjects

  x_min <- min(x, na.rm = TRUE)
  mu0 <- stats::median(y[x == x_min], na.rm = TRUE)
  # Guard mu0 into the open interval so the inversion is finite.
  mu0 <- min(max(mu0, 1e-3), 1 - 1e-3)

  if (x_min <= 0) {
    # Cannot invert discounting function at x=0 (division by zero); use safe default.
    k0 <- 0.01
  } else if (identical(equation, "exponential")) {
    # mu = exp(-k * x_min)  =>  k = -log(mu) / x_min
    k0 <- -log(mu0) / x_min
  } else {
    # mazur: mu = 1 / (1 + k * x_min)  =>  k = (1/mu - 1) / x_min
    k0 <- (1 / mu0 - 1) / x_min
  }
  if (!is.finite(k0) || k0 <= 0) k0 <- 0.01
  log_k0 <- log(k0)

  beta_k <- rep(0, p)
  beta_k[1] <- log_k0

  log_aux <- if (identical(family, "gaussian")) log(0.1) else log(8)

  list(
    beta_k = beta_k,
    log_sigma_u = log(0.5),
    log_aux = log_aux,
    u = matrix(0, nrow = n_subjects, ncol = 1L)
  )
}


# ==============================================================================
# P.4  Optimizer dispatch + bounds expansion (verbatim port from beezdemand)
# ==============================================================================

#' Expand partial optimizer bounds to the full parameter vector
#'
#' Ported verbatim from beezdemand `.expand_bounds`. Given a named vector of
#' user-specified bounds (possibly partial, possibly `NULL`), expands it to the
#' full length of `par_names` by filling unspecified positions with
#' `default_val`. Repeated parameter names (e.g. `beta_k`, `beta_k`) are each
#' filled when the corresponding name appears in `bounds`.
#'
#' @param bounds Named numeric vector of user-specified bounds, or `NULL`.
#' @param par_names Character vector of optimizer parameter names (from
#'   `names(obj$par)`). May contain repeated names for vector parameters.
#' @param default_val Default bound value: `-Inf` for lower, `Inf` for upper.
#' @return Numeric vector of length `length(par_names)`.
#' @keywords internal
.expand_bounds <- function(bounds, par_names, default_val) {
  if (is.null(bounds)) return(rep(default_val, length(par_names)))
  result <- rep(default_val, length(par_names))
  names(result) <- par_names
  for (nm in names(bounds)) {
    idx <- which(par_names == nm)
    if (length(idx) > 0) {
      result[idx] <- bounds[nm]
    } else {
      warning("Bounds specified for unknown parameter '", nm, "' (ignored)",
              call. = FALSE)
    }
  }
  result
}


#' Run a single TMB optimization (nlminb or L-BFGS-B), family-agnostic
#'
#' Ported verbatim from beezdemand `.tmb_run_optimizer` (renamed only).
#' Dispatches to `nlminb` (default) or `optim(method = "L-BFGS-B")` and
#' normalizes the return value so downstream code always sees identical field
#' names regardless of optimizer: `$par`, `$objective`, `$convergence`,
#' `$message` (guaranteed non-`NULL` character).
#'
#' Optimizer warnings are captured via `withCallingHandlers` and returned in
#' `$warnings`; errors yield `convergence = 99L, objective = Inf`.
#'
#' @param obj TMB objective function object (from `TMB::MakeADFun`), with
#'   `$fn`, `$gr`, and `$par`.
#' @param start Named numeric vector of starting parameter values.
#' @param tmb_control Merged control list; must have fields `optimizer`,
#'   `iter_max`, `eval_max`, `rel_tol`, `lower`, `upper`, `trace`.
#' @param user_specified Character vector of fields the user explicitly set in
#'   `tmb_control` (governs `trace` precedence).
#' @param verbose Integer verbosity level (0 = silent, 1 = progress,
#'   2 = full trace).
#' @return A list with `opt` (list of `par`, `objective`, `convergence`,
#'   `message`) and `warnings` (character vector).
#' @keywords internal
.dd_tmb_run_optimizer <- function(obj, start, tmb_control, user_specified,
                                   verbose) {
  optimizer <- tmb_control$optimizer
  iter_max  <- tmb_control$iter_max
  eval_max  <- tmb_control$eval_max
  rel_tol   <- tmb_control$rel_tol

  # Effective trace: user-specified takes precedence, then verbose.
  trace <- if ("trace" %in% user_specified) {
    as.integer(tmb_control$trace)
  } else if (verbose >= 2) {
    1L
  } else {
    0L
  }

  # Expand bounds to the full parameter length.
  par_names <- names(start)
  lower <- .expand_bounds(tmb_control$lower, par_names, -Inf)
  upper <- .expand_bounds(tmb_control$upper, par_names,  Inf)

  opt_warnings <- character(0)

  if (optimizer == "nlminb") {
    opt <- tryCatch(
      withCallingHandlers(
        stats::nlminb(
          start      = start,
          objective  = obj$fn,
          gradient   = obj$gr,
          lower      = lower,
          upper      = upper,
          control    = list(
            eval.max = eval_max,
            iter.max = iter_max,
            rel.tol  = rel_tol,
            trace    = trace
          )
        ),
        warning = function(w) {
          opt_warnings <<- c(opt_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        list(
          par         = start,
          objective   = Inf,
          convergence = 99L,
          message     = conditionMessage(e)
        )
      }
    )
  } else {
    # L-BFGS-B via stats::optim
    opt <- tryCatch(
      {
        raw <- withCallingHandlers(
          stats::optim(
            par     = start,
            fn      = obj$fn,
            gr      = obj$gr,
            method  = "L-BFGS-B",
            lower   = lower,
            upper   = upper,
            control = list(
              maxit = iter_max,
              trace = trace
            )
          ),
          warning = function(w) {
            opt_warnings <<- c(opt_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )
        # Normalize optim's $value/$convergence to nlminb field names.
        list(
          par         = raw$par,
          objective   = raw$value,
          convergence = raw$convergence,
          message     = raw$message %||% "maximum iterations reached"
        )
      },
      error = function(e) {
        list(
          par         = start,
          objective   = Inf,
          convergence = 99L,
          message     = conditionMessage(e)
        )
      }
    )
  }

  # Guarantee $message is a non-NULL character string.
  if (is.null(opt$message)) opt$message <- "unknown"

  list(opt = opt, warnings = opt_warnings)
}


# ==============================================================================
# P.5  Multi-start + degenerate-optimum guard
# ==============================================================================

# phi floor for the sltb family: imposed as an OPTIMIZER LOWER BOUND on
# log_aux (log_aux >= log(.dd_phi_min)) so the optimizer can never reach the
# phi->0 / k->inf collapse (spec section 4.8). Overridable via
# tmb_control$lower$log_aux; the user value then wins.
.dd_phi_min <- 0.1

# Sanity NaN/blowup guard on the log-k intercept ONLY (not a phi rejection).
# A fit where |beta_k[1]| > 20 implies k outside [exp(-20), exp(20)] ~ [2e-9, 5e8],
# which is numerically nonsensical for any discounting experiment.
.dd_beta_k_abs_max <- 20


#' Multi-start optimization for the TMB discounting model
#'
#' Builds 3 starting sets (data-driven, low-k, high-k) and runs each through
#' [.dd_tmb_run_optimizer()]. Keeps the lowest **finite** nll among candidates
#' passing a NaN/blowup sanity check on `beta_k[1]`.
#'
#' **Degenerate-optimum guard (phi->0 / k->inf, spec section 4.8):**
#' For `family == "sltb"` (i.e. `tmb_data$family == 0L`) a **lower bound on
#' `log_aux`** of `log(.dd_phi_min)` (= `log(0.1)`) is inserted into
#' `tmb_control$lower` before every optimizer call so the optimizer can never
#' walk into the degenerate phi->0 optimum. The bound is overridable: if the
#' caller already provides a `log_aux` entry in `tmb_control$lower` that value
#' wins (merged before the floor check). Genuine low-precision fits (true phi
#' above the floor, e.g. phi = 2) are retained — the bound only blocks the
#' pathological phi->0 sink. There is **no post-hoc phi-based rejection**.
#'
#' If every candidate trips the blowup guard the lowest-nll fit overall is
#' returned with a warning. If every candidate fails entirely (error or
#' non-finite nll) the function stops.
#'
#' @param tmb_data TMB data list from [.dd_tmb_build_tmb_data()].
#' @param start_values Default starting list from [.dd_tmb_default_starts()].
#' @param tmb_control Merged control list (may carry user `lower`/`upper`).
#' @param user_specified Character vector of user-set `tmb_control` fields.
#' @param verbose Integer verbosity level.
#' @return A list with elements `obj`, `opt`, `nll`, `start_idx`,
#'   `opt_warnings`.
#' @keywords internal
.dd_tmb_multi_start <- function(tmb_data, start_values,
                                 tmb_control, user_specified, verbose) {
  is_sltb <- identical(tmb_data$family, 0L)

  # Impose the phi floor as a log_aux lower bound for sltb, unless the user
  # already supplied a log_aux lower bound (their value wins).
  if (is_sltb) {
    user_lower        <- tmb_control$lower
    has_user_log_aux  <- !is.null(user_lower) &&
                           "log_aux" %in% names(user_lower)
    if (!has_user_log_aux) {
      floor_bound          <- c(log_aux = log(.dd_phi_min))
      tmb_control$lower    <- c(user_lower, floor_bound)
    }
  }

  # 3 starting sets: data-driven (default), low-k, high-k.
  start_sets      <- vector("list", 3L)
  start_sets[[1]] <- start_values                    # data-driven

  s2               <- start_values                   # low-k
  s2$beta_k[1]    <- start_values$beta_k[1] - 1.5
  s2$log_sigma_u  <- log(0.3)
  start_sets[[2]] <- s2

  s3               <- start_values                   # high-k
  s3$beta_k[1]    <- start_values$beta_k[1] + 1.5
  s3$log_sigma_u  <- log(0.8)
  start_sets[[3]] <- s3

  # Sanity NaN/blowup guard on the log-k intercept ONLY (no phi rejection: the
  # log_aux lower bound already prevents the degenerate phi->0 optimum).
  .is_blowup <- function(opt) {
    par   <- opt$par
    beta0 <- par[names(par) == "beta_k"][1]
    !is.finite(beta0) || abs(beta0) > .dd_beta_k_abs_max
  }

  best_kept_nll <- Inf
  best_kept     <- NULL
  best_any_nll  <- Inf
  best_any      <- NULL

  for (s in seq_along(start_sets)) {
    starts_i <- start_sets[[s]]

    result <- tryCatch({
      obj_i <- TMB::MakeADFun(
        data       = tmb_data,
        parameters = starts_i,
        random     = "u",
        DLL        = "beezdiscounting",
        silent     = verbose < 2
      )
      opt_res_i <- .dd_tmb_run_optimizer(
        obj_i, obj_i$par, tmb_control, user_specified, verbose
      )
      list(
        obj         = obj_i,
        opt         = opt_res_i$opt,
        nll         = opt_res_i$opt$objective,
        start_idx   = s,
        opt_warnings = opt_res_i$warnings
      )
    }, error = function(e) {
      if (verbose >= 2) {
        message(sprintf("  Start set %d failed: %s", s, e$message))
      }
      NULL
    })

    if (is.null(result) || !is.finite(result$nll)) next

    if (result$nll < best_any_nll) {
      best_any_nll <- result$nll
      best_any     <- result
    }
    if (!.is_blowup(result$opt) && result$nll < best_kept_nll) {
      best_kept_nll <- result$nll
      best_kept     <- result
    }
  }

  if (is.null(best_any)) {
    stop(
      "All starting value sets failed. ",
      "Check data quality or try different start values.",
      call. = FALSE
    )
  }

  if (!is.null(best_kept)) {
    best_result <- best_kept
  } else {
    best_result <- best_any
    warning(
      "All multi-start fits tripped the beta_k sanity guard ",
      "(k -> Inf / non-finite intercept); returning the lowest-nll fit. ",
      "Inspect data for boundary-heavy subjects.",
      call. = FALSE
    )
  }

  if (verbose >= 1) {
    message(sprintf(
      "  Multi-start: best NLL = %.2f (start set %d of %d)",
      best_result$nll, best_result$start_idx, length(start_sets)
    ))
  }

  best_result
}


#' Extract estimates from a TMB discounting fit
#'
#' Computes `TMB::sdreport`, gates on `isTRUE(sdr$pdHess)` (warn, never abort),
#' and renames the generic auxiliary scalar `log_aux` to `log_phi` (sltb) or
#' `log_sigma_e` (gaussian) in the returned `coefficients`/`se`.
#'
#' @param obj TMB objective object.
#' @param opt Normalized optimizer result (from [.dd_tmb_run_optimizer()]).
#' @param n_subjects Integer number of subjects.
#' @param family One of "sltb", "gaussian".
#' @param verbose Integer verbosity.
#' @return list(coefficients, se, sdr, variance_components, u_hat, hessian_pd).
#' @keywords internal
.dd_tmb_extract_estimates <- function(obj, opt, n_subjects, family, verbose = 1) {
  sdr <- tryCatch(
    TMB::sdreport(obj),
    error = function(e1) {
      sdr2 <- tryCatch(
        TMB::sdreport(obj, getJointPrecision = FALSE),
        error = function(e2) NULL
      )
      if (is.null(sdr2) && verbose >= 1) {
        warning("Standard error computation failed: ", e1$message)
      }
      sdr2
    }
  )

  hessian_pd <- NA
  if (!is.null(sdr)) {
    hessian_pd <- isTRUE(sdr$pdHess)
    if (!hessian_pd && verbose >= 1) {
      cli::cli_warn(c(
        "!" = "Hessian is not positive definite ({.code pdHess = FALSE}).",
        "i" = "Standard errors, p-values, and confidence intervals may be unreliable.",
        "i" = "Consider simplifying the model or checking data quality."
      ))
    }
  }

  par_full <- opt$par
  par_names <- names(par_full)

  coefficients <- par_full
  se_vec <- rep(NA_real_, length(par_full))
  names(se_vec) <- par_names

  if (!is.null(sdr)) {
    fixed_summary <- summary(sdr, "fixed")
    .fill_vector_se <- function(name) {
      idx <- which(par_names == name)
      if (length(idx) == 0L) return(invisible(NULL))
      rows <- fixed_summary[rownames(fixed_summary) == name, , drop = FALSE]
      if (nrow(rows) == length(idx)) se_vec[idx] <<- rows[, "Std. Error"]
    }
    .fill_vector_se("beta_k")
    .fill_vector_se("log_sigma_u")
    .fill_vector_se("log_aux")

    re_summary <- tryCatch(summary(sdr, "random"), error = function(e) NULL)
    if (!is.null(re_summary)) {
      u_hat <- matrix(re_summary[, "Estimate"], nrow = n_subjects, ncol = 1L)
    } else {
      u_hat <- matrix(0, nrow = n_subjects, ncol = 1L)
    }
  } else {
    u_hat <- matrix(0, nrow = n_subjects, ncol = 1L)
  }

  variance_components <- NULL
  if (!is.null(sdr)) {
    variance_components <- tryCatch(summary(sdr, "report"), error = function(e) NULL)
  }

  # Rename the generic auxiliary scalar in both coefficients and se.
  aux_name <- if (identical(family, "gaussian")) "log_sigma_e" else "log_phi"
  names(coefficients)[names(coefficients) == "log_aux"] <- aux_name
  names(se_vec)[names(se_vec) == "log_aux"] <- aux_name

  list(
    coefficients = coefficients,
    se = se_vec,
    sdr = sdr,
    variance_components = variance_components,
    u_hat = u_hat,
    hessian_pd = hessian_pd
  )
}


#' Compute subject-specific discounting parameters
#'
#' Reconstructs each subject's `k_i = exp(beta_k[1] + sigma_u * u_i)` using the
#' non-centered predictor that matches the C++ template
#' (`log_k_i = X.row(i)*beta_k + sigma_u * u(subj,0)`). For the intercept-only
#' MVP, `Xbeta` equals the intercept for all subjects; subject k differs via
#' `u_i`. The auxiliary scalar `phi` is population-level in the MVP, so it is
#' **not** a subject-level parameter and is never returned here (for either
#' family).
#'
#' @param coefficients Named coefficient vector (with `beta_k`, `log_sigma_u`,
#'   and `log_phi` or `log_sigma_e`).
#' @param u_hat Matrix `n_subjects` x 1 of standardized random effects.
#' @param subject_levels Character vector of subject ids (length n_subjects).
#' @param equation One of "mazur", "exponential" (reserved; k is equation-free).
#' @param family One of "sltb", "gaussian".
#' @return data.frame(id, u_i, k) — no phi column.
#' @note For factor designs with multiple `beta_k` columns, only `beta_k[1]`
#'   (the population intercept) is used when computing per-subject k. Between-
#'   subject factor contributions to the log-k linear predictor are therefore
#'   ignored. This is the correct MVP behavior for the intercept-only design
#'   (single `beta_k`), but will under-estimate k for subjects in non-reference
#'   factor groups in multi-factor designs. Use `predict()` for cell-level values
#'   once factor support is added in a future phase.
#' @keywords internal
.dd_tmb_compute_subject_pars <- function(coefficients, u_hat, subject_levels,
                                         equation, family) {
  beta_k <- unname(coefficients[names(coefficients) == "beta_k"])
  beta0 <- beta_k[1]
  sigma_u <- exp(unname(coefficients[["log_sigma_u"]]))
  u_i <- as.numeric(u_hat[, 1L])

  log_k_i <- beta0 + sigma_u * u_i
  k_i <- exp(log_k_i)

  data.frame(
    id = subject_levels,
    u_i = u_i,
    k = k_i,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# P.8  Public entry point: fit_dd_tmb()
# ==============================================================================

#' Fit an indifference-point mixed-effects discounting model via TMB
#'
#' Fits a 1-parameter discounting model (Mazur hyperbolic or exponential) with a
#' random intercept on `log k`, between-subject fixed effects, and either an
#' SLT-beta or Gaussian observation family, using Template Model Builder for
#' exact AD + Laplace approximation.
#'
#' @param data Long data frame with subject id, delay, and indifference
#'   proportion columns.
#' @param y_var,x_var,id_var Column names (defaults `"y"`, `"x"`, `"id"`).
#' @param equation One of `"mazur"`, `"exponential"`.
#' @param family Observation family: `"sltb"` (default) or `"gaussian"`.
#' @param random_effects RE formula (MVP: `k ~ 1`).
#' @param factors Character vector of between-subject factor names.
#' @param factor_interaction Logical; include a pairwise factor interaction.
#' @param continuous_covariates Character vector of covariate names.
#' @param ll Optional larger-later reward for `amount`-scale coercion.
#' @param response_scale One of `"proportion"`, `"percent"`, `"amount"`.
#' @param start_values Optional named list overriding defaults.
#' @param tmb_control Optimizer control list.
#' @param multi_start Logical; if `TRUE` (default), run the 3-set guarded
#'   multi-start.
#' @param verbose Integer verbosity (0 silent, 1 progress, 2 debug).
#' @param ... Reserved.
#' @return An object of class `beezdiscounting_tmb` with components:
#'   \describe{
#'     \item{call}{The matched call.}
#'     \item{opt}{Normalized optimizer result (`par`, `objective`,
#'       `convergence`, `message`).}
#'     \item{model}{List of `coefficients`, `se`, and `variance_components`.}
#'     \item{sdr}{TMB `sdreport` object (or `NULL` if SE computation failed).}
#'     \item{hessian_pd}{Logical positive-definiteness of the Hessian.}
#'     \item{param_info}{Model metadata (equation, family, dimensions, factor
#'       spec, parsed random effects).}
#'     \item{formula_details}{Fixed-effect design (`X`, `rhs`, `contrasts`).}
#'     \item{subject_pars}{Data frame of subject-level `id`, `u_i`, `k`.}
#'     \item{loglik, AIC, BIC}{Fit statistics.}
#'     \item{converged, se_available}{Convergence / SE-availability flags.}
#'     \item{opt_warnings}{Character vector of optimizer warnings.}
#'     \item{data}{The single filtered model frame (id/x/y + retained design
#'       columns), row-aligned with the design matrix.}
#'     \item{data_all}{The validated frame before complete-casing.}
#'     \item{coercion_info}{Scale-coercion/clamping audit list.}
#'   }
#' @examples
#' \donttest{
#' # Small two-subject long-format indifference-point data frame.
#' dd <- data.frame(
#'   id = rep(c("s1", "s2"), each = 5),
#'   x  = rep(c(7, 30, 180, 365, 730), times = 2),
#'   y  = c(0.95, 0.80, 0.45, 0.30, 0.15,
#'          0.90, 0.70, 0.40, 0.25, 0.10)
#' )
#' fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb", verbose = 0)
#' exp(fit$model$coefficients[["beta_k"]])  # population k
#' }
#' @export
fit_dd_tmb <- function(data,
                       y_var = "y", x_var = "x", id_var = "id",
                       equation = c("mazur", "exponential"),
                       family = c("sltb", "gaussian"),
                       random_effects = k ~ 1,
                       factors = NULL,
                       factor_interaction = FALSE,
                       continuous_covariates = NULL,
                       ll = NULL,
                       response_scale = c("proportion", "percent", "amount"),
                       start_values = NULL,
                       tmb_control = list(iter_max = 1000, eval_max = 2000),
                       multi_start = TRUE,
                       verbose = 1,
                       ...) {
  cl <- match.call()
  equation <- match.arg(equation)
  family <- match.arg(family)
  response_scale <- match.arg(response_scale)

  # The union of factor + covariate columns is the set of EXTRA modeling
  # columns that must travel with id/x/y through one complete-case pass
  # (row-coherence). The validator retains them, prepare_data complete-cases
  # them ONCE, and the design + TMB arrays are derived from that single frame.
  extra_cols <- unique(c(factors, continuous_covariates))

  # 1. Validate + coerce/clamp (warns; errors on missing/NA y). Retains the
  #    factor/covariate columns alongside canonical id/x/y.
  validated <- .dd_validate_ip(data, y_var = y_var, x_var = x_var,
                               id_var = id_var, ll = ll,
                               extra_cols = extra_cols,
                               response_scale = response_scale)
  long <- validated$data            # canonical id/x/y + retained extras
  coercion_info <- validated$coercion_info

  # 2. Random-effects normalization (MVP: single intercept-only block).
  re_norm <- .dd_normalize_re(random_effects, data = long)
  n_random_effects <- 1L

  # 3. Prepare data: ONE complete-case pass over id/x/y + extra_cols, building
  #    the 0-indexed subject_id. prepared$data is the single filtered model
  #    frame (id/x/y + retained extras) that everything else derives from.
  prepared <- .dd_tmb_prepare_data(long, y_var = "y", x_var = "x",
                                   id_var = "id", extra_cols = extra_cols)

  # 4. Fixed-effect design for log k, built on the SAME filtered frame
  #    (prepared$data) so X's rows align 1:1 with the prepared y/x/subject_id.
  design <- .dd_tmb_build_design(prepared$data, factors = factors,
                                 factor_interaction = factor_interaction,
                                 continuous_covariates = continuous_covariates)

  # 5. TMB data + default starts.
  tmb_data <- .dd_tmb_build_tmb_data(prepared, design, equation, family)
  default_starts <- .dd_tmb_default_starts(prepared, design, family, equation)
  if (!is.null(start_values)) {
    for (nm in names(start_values)) {
      if (nm %in% names(default_starts)) default_starts[[nm]] <- start_values[[nm]]
    }
  }

  # 6. Merge control defaults.
  default_control <- list(
    iter_max = 1000, eval_max = 2000, optimizer = "nlminb",
    rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0
  )
  user_specified <- names(tmb_control)
  tmb_control <- utils::modifyList(default_control, tmb_control)

  if (verbose >= 1) {
    message(sprintf("Fitting TMB mixed-effects discounting model (%s, %s)...",
                    equation, family))
    message(sprintf("  Subjects: %d, Observations: %d",
                    prepared$n_subjects, prepared$n_obs))
  }

  # 7. Optimize (multi-start or single).
  opt_warnings <- character(0)
  if (isTRUE(multi_start)) {
    result <- .dd_tmb_multi_start(tmb_data, default_starts, tmb_control,
                                  user_specified, verbose)
    obj <- result$obj
    opt <- result$opt
    opt_warnings <- result$opt_warnings %||% character(0)
  } else {
    obj <- TMB::MakeADFun(tmb_data, default_starts, random = "u",
                          DLL = "beezdiscounting", silent = verbose < 2)
    opt_res <- .dd_tmb_run_optimizer(obj, obj$par, tmb_control,
                                     user_specified, verbose)
    opt <- opt_res$opt
    opt_warnings <- opt_res$warnings
  }

  converged <- isTRUE(opt$convergence == 0)
  try(obj$fn(opt$par), silent = TRUE)

  # 8. Extract estimates (sdreport, pdHess gate, log_aux rename).
  estimates <- .dd_tmb_extract_estimates(obj, opt,
                                         n_subjects = prepared$n_subjects,
                                         family = family, verbose = verbose)

  # 9. Subject-specific parameters (id/u_i/k; no phi -- phi is population-level).
  subject_pars <- .dd_tmb_compute_subject_pars(
    coefficients = estimates$coefficients,
    u_hat = estimates$u_hat,
    subject_levels = prepared$subject_levels,
    equation = equation,
    family = family
  )

  # 10. Likelihood / IC.
  nll <- opt$objective
  loglik <- -nll
  n_fixed_params <- length(opt$par)
  aic <- 2 * nll + 2 * n_fixed_params
  bic <- 2 * nll + n_fixed_params * log(prepared$n_obs)

  has_phi <- identical(family, "sltb")

  result_obj <- structure(
    list(
      call = cl,
      opt = opt,
      model = list(
        coefficients = estimates$coefficients,
        se = estimates$se,
        variance_components = estimates$variance_components
      ),
      sdr = estimates$sdr,
      hessian_pd = estimates$hessian_pd,
      param_info = list(
        equation = equation,
        family = family,
        has_phi = has_phi,
        n_obs = prepared$n_obs,
        n_subjects = prepared$n_subjects,
        n_random_effects = n_random_effects,
        subject_levels = prepared$subject_levels,
        id_var = id_var,
        x_var = x_var,
        y_var = y_var,
        # Factor metadata persisted so predict()/emmeans can rebuild a
        # column-aligned design via build_fixed_rhs + stored contrasts (B3).
        factors = factors,
        factor_interaction = factor_interaction,
        continuous_covariates = continuous_covariates,
        random_effects_parsed = re_norm
      ),
      # Store rhs + contrasts so .dd_build_emm_ref_grid (E.1) and predict (M.3)
      # rebuild the design through the SAME build_fixed_rhs route as the fit.
      formula_details = list(X = design$X, rhs = design$rhs,
                             contrasts = design$contrasts),
      subject_pars = subject_pars,
      loglik = loglik,
      AIC = aic,
      BIC = bic,
      converged = converged,
      se_available = !is.null(estimates$sdr),
      opt_warnings = opt_warnings,
      # fit$data is the SINGLE filtered model frame (id/x/y + retained factor/
      # covariate columns), row-aligned with X and the prepared arrays (B3).
      data = prepared$data,
      data_all = long,
      coercion_info = coercion_info
    ),
    class = "beezdiscounting_tmb"
  )

  if (verbose >= 1) {
    if (converged) {
      message(sprintf("  Converged (NLL = %.2f). Done.", nll))
    } else {
      message(sprintf("  WARNING: Did not converge (code %s: %s).",
                      opt$convergence, opt$message))
    }
  }

  result_obj
}
