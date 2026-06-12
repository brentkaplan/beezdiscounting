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

  # B3: reject a rank-deficient design. An empty factor cell (e.g. an absent
  # interaction combination in an unbalanced design) yields an all-zero, aliased
  # column; TMB would then "estimate" an unidentified coefficient and the EMM
  # marginalization would average over a non-estimable cell -- a confident but
  # meaningless result. qr()'s default tolerance (1e-7) matches base `lm`'s
  # aliasing detection; the pivot names the dependent columns.
  if (ncol(X) > 0L) {
    qr_x <- qr(X)
    if (qr_x$rank < ncol(X)) {
      aliased <- colnames(X)[qr_x$pivot[(qr_x$rank + 1L):ncol(X)]]
      cli::cli_abort(c(
        "The fixed-effect design for {.code log k} is rank-deficient \\
         ({qr_x$rank}/{ncol(X)} columns independent).",
        "x" = "Non-estimable (aliased) column{?s}: {.val {aliased}}.",
        "i" = "This usually means an empty factor cell or a redundant \\
               interaction. Drop the interaction ({.code factor_interaction = \\
               FALSE}), collapse empty cells, or remove the collinear predictor."
      ))
    }
  }

  list(X = X, rhs = rhs, contrasts = attr(X, "contrasts"))
}


#' Build the TMB data list for MixedDiscounting
#'
#' @param prepared Output from [.dd_tmb_prepare_data()].
#' @param design Output from [.dd_tmb_build_design()].
#' @param equation One of "mazur", "exponential".
#' @param family One of "sltb", "gaussian".
#' @param n_re Number of random-effect intercepts: `1L` (on `log k`) or `2L`
#'   (joint `(log k, log phi)`). Defaults to `1L`.
#' @return A list whose names match the C++ `DATA_*` macros: `model`, `y`, `x`,
#'   `subject_id` (0-indexed integer), `X`, `eqn_type`, `family`, `n_obs`,
#'   `n_subjects`, `n_re`.
#' @keywords internal
.dd_tmb_build_tmb_data <- function(prepared, design, equation, family,
                                   n_re = 1L, re2_target = 0L) {
  eqn_type <- switch(equation,
    mazur = 0L,
    exponential = 1L,
    `green-myerson` = 2L,
    rachlin = 3L,
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
    n_subjects = as.integer(prepared$n_subjects),
    n_re = as.integer(n_re),
    re2_target = as.integer(re2_target)
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
#'   `log_aux`, `log_s` (always present; held fixed for 1-parameter equations),
#'   `u` (matrix `n_subjects` x 1), plus the joint 2-RE blocks `log_sd_re`
#'   (length 2), `cor_re` (length 1), and `b` (matrix `n_subjects` x 2). The
#'   shared C++ template declares every block; the map fixes the inactive ones
#'   per `n_re`.
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
    log_s = 0,                       # s = 1 start; map-fixed for 1-param eqns
    u = matrix(0, nrow = n_subjects, ncol = 1L),
    # Joint 2-RE (log k, log phi) blocks; always emitted so the shared template
    # has every parameter. The map fixes these for a 1-RE fit (n_re == 1L).
    log_sd_re = c(log(0.5), log(0.5)),
    cor_re = 0,
    b = matrix(0, nrow = n_subjects, ncol = 2L)
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

# Sanity NaN/blowup guard. A fitted |log k| > 20 implies k outside
# [exp(-20), exp(20)] ~ [2e-9, 5e8], numerically nonsensical for any discounting
# experiment.
.dd_beta_k_abs_max <- 20


# B7: blow-up predicate on the FULL fitted log-k linear predictor
# (eta_i = X.row(i) %*% beta_k), not just the intercept -- a degenerate
# non-reference level / interaction / covariate slope can push some subject's
# log k past the bound while beta_k[1] stays sane. Returns TRUE on a beta/X shape
# mismatch, non-finite beta or eta, or max |eta| over the bound.
.dd_logk_blowup <- function(opt, design_X, max_abs = .dd_beta_k_abs_max) {
  par  <- opt$par
  beta <- unname(par[names(par) == "beta_k"])
  x_mat <- as.matrix(design_X)
  if (length(beta) != ncol(x_mat)) return(TRUE)
  if (any(!is.finite(beta))) return(TRUE)
  eta <- as.numeric(x_mat %*% beta)
  if (any(!is.finite(eta))) return(TRUE)
  max(abs(eta)) > max_abs
}


#' Impose the SLT-beta phi floor as a log_aux lower bound
#'
#' For `family == "sltb"` the degenerate phi->0 / k->inf optimum (spec section
#' 4.8) is blocked by a **lower bound on `log_aux`** of `log(.dd_phi_min)`
#' inserted into `tmb_control$lower`. This must apply to BOTH optimizer paths
#' (single fit and multi-start), so the logic lives in one shared helper. The
#' bound is overridable: if the caller already supplies a `log_aux` entry in
#' `tmb_control$lower`, their value wins (the floor is not added).
#'
#' @param tmb_control Merged control list (may carry user `lower`/`upper`).
#' @param family One of "sltb", "gaussian".
#' @return The (possibly augmented) `tmb_control` list.
#' @keywords internal
.dd_apply_phi_floor <- function(tmb_control, family) {
  if (!identical(family, "sltb")) return(tmb_control)
  user_lower       <- tmb_control$lower
  has_user_log_aux <- !is.null(user_lower) && "log_aux" %in% names(user_lower)
  if (!has_user_log_aux) {
    tmb_control$lower <- c(user_lower, log_aux = log(.dd_phi_min))
  }
  tmb_control
}


# Soft (wide) optimizer bounds on log_s for the 2-parameter equations; a
# numerical guard, NOT a constraint of scientific interest (no hard s <= 1 cap).
.dd_s_log_lower <- log(0.05)
.dd_s_log_upper <- log(20)


#' Build the TMB `map` for the discounting shape + random-effect blocks
#'
#' Generalizes the old `has_s`-only map to cover the joint 2-RE blocks. The
#' shared C++ template declares EVERY parameter (`log_s`, `log_sigma_u`/`u`,
#' `log_sd_re`/`cor_re`/`b`); this map fixes the inactive ones so a given fit
#' only estimates the blocks it should:
#'
#' - `log_s` is held fixed (`factor(NA)`) for the 1-parameter equations
#'   (`!has_s`) and freed for the 2-parameter equations (`has_s`).
#' - `n_re == 1L` (single intercept on `log k`): fix `log_sd_re`, `cor_re` (and
#'   `b` is fixed by [.dd_tmb_finalize_map()], which needs the starts to size
#'   the factor).
#' - `n_re == 2L` (joint `(log k, log phi)`): fix `log_sigma_u` (and `u` via
#'   [.dd_tmb_finalize_map()]); fix `cor_re` iff `covariance == "pdDiag"`.
#'
#' The returned value is finalized by [.dd_tmb_finalize_map()] and passed
#' verbatim to EVERY `TMB::MakeADFun()` call so a fit can never estimate an
#' unidentified block (singular Hessian, wrong df).
#'
#' @param has_s Logical; `TRUE` for green-myerson / rachlin.
#' @param n_re Number of random-effect intercepts (`1L` or `2L`).
#' @param covariance `"pdSymm"` (cor_re free) or `"pdDiag"` (cor_re fixed); only
#'   consulted when `n_re == 2L`.
#' @return A (possibly empty) named list of `factor(NA)` map entries.
#' @keywords internal
.dd_tmb_build_map <- function(has_s, n_re = 1L, covariance = "pdSymm") {
  map <- list()
  if (!isTRUE(has_s)) map$log_s <- factor(NA)
  if (n_re == 1L) {
    map$log_sd_re <- factor(c(NA, NA))
    map$cor_re <- factor(NA)
    # The inactive matrix (b for 1-RE; u for 2-RE) is fixed at 0 by
    # .dd_tmb_finalize_map(), which needs the starts to size the factor().
  } else {
    map$log_sigma_u <- factor(NA)
    if (identical(covariance, "pdDiag")) map$cor_re <- factor(NA)
  }
  map
}


#' Fix the inactive random-effect matrix in the TMB map
#'
#' `u` (1-RE) and `b` (2-RE) are alternately the integrated random block vs a
#' fixed-at-zero block. This appends the `factor(NA)` fix for whichever matrix
#' is inactive, sized from its starts. Applied right before each
#' `TMB::MakeADFun()`.
#'
#' @param map The map from [.dd_tmb_build_map()].
#' @param starts The starting-value list (for the matrix lengths).
#' @param n_re Number of random-effect intercepts (`1L` or `2L`).
#' @return The augmented map.
#' @keywords internal
.dd_tmb_finalize_map <- function(map, starts, n_re) {
  if (n_re == 1L) {
    map$b <- factor(rep(NA, length(starts$b)))
  } else {
    map$u <- factor(rep(NA, length(starts$u)))
  }
  map
}


#' Assert the exact free-parameter block set at the map seam
#'
#' Pins the free-block set for each `(n_re, covariance, has_s)` combination so a
#' leaked block (an unidentified parameter slipping into the optimizer) or a
#' missing block (a free coefficient accidentally mapped out) is caught
#' immediately rather than surfacing as a silent fit pathology.
#'
#' @param par_names Character vector of optimizer parameter names
#'   (`names(opt$par)`).
#' @param expected Character vector of expected free-block names.
#' @param context Short label for the error message.
#' @return Invisibly `TRUE`; stops on mismatch.
#' @keywords internal
.dd_assert_free_params <- function(par_names, expected, context) {
  got <- sort(unique(par_names))
  exp <- sort(unique(expected))
  if (!identical(got, exp)) {
    stop("internal: free-parameter blocks for ", context, " do not match the ",
         "expected map (got {", paste(got, collapse = ", "), "}; expected {",
         paste(exp, collapse = ", "), "}); check map threading.", call. = FALSE)
  }
  invisible(TRUE)
}


# 2-RE degeneracy guard (spec section 4): a wide lower+upper bound on log_sd_re
# so the RE variance cannot run away and drive a boundary subject's phi_i -> 0
# (the log_aux floor only guards the POPULATION phi). Mirrors .dd_apply_s_bounds:
# only adds bounds when n_re == 2 and the caller has not already set them.
.dd_re_log_sd_lower <- log(1e-3)
.dd_re_log_sd_upper <- log(5)

#' Impose wide optimizer bounds on log_sd_re for a 2-RE fit
#'
#' Adds lower/upper `log_sd_re` bounds (`[log(1e-3), log(5)]`) only when
#' `n_re == 2L` and the caller has not already set a `log_sd_re` bound. Each
#' bound is duplicated (the vector has length 2). Keeps Sigma positive-definite
#' and stops an absurd RE-variance estimate from collapsing a boundary subject's
#' `phi_i`.
#'
#' @param tmb_control Merged control list (may carry user `lower`/`upper`).
#' @param n_re Number of random-effect intercepts (`1L` or `2L`).
#' @return The (possibly augmented) `tmb_control` list.
#' @keywords internal
.dd_apply_re_bounds <- function(tmb_control, n_re) {
  if (!identical(as.integer(n_re), 2L)) return(tmb_control)
  lo <- tmb_control$lower
  if (is.null(lo) || !("log_sd_re" %in% names(lo))) {
    # Both log_sd_re entries must be equal: .expand_bounds recycles only the
    # first matched value across all log_sd_re optimizer positions, so differing
    # per-axis values would be silently dropped.
    tmb_control$lower <- c(lo, log_sd_re = .dd_re_log_sd_lower,
                           log_sd_re = .dd_re_log_sd_lower)
  }
  up <- tmb_control$upper
  if (is.null(up) || !("log_sd_re" %in% names(up))) {
    tmb_control$upper <- c(up, log_sd_re = .dd_re_log_sd_upper,
                           log_sd_re = .dd_re_log_sd_upper)
  }
  tmb_control
}


#' Impose wide optimizer bounds on log_s for the 2-parameter equations
#'
#' Mirrors [.dd_apply_phi_floor()]: adds `log_s` lower/upper entries to
#' `tmb_control$lower`/`$upper` (wide numerical guards `[log(0.05), log(20)]`)
#' only when `has_s` and only when the caller has not already set a `log_s`
#' bound (user value wins). For 1-parameter equations `log_s` is mapped out of
#' the optimizer vector, so no bound is added (and `.expand_bounds` would
#' otherwise warn about an unknown parameter).
#'
#' @param tmb_control Merged control list (may carry user `lower`/`upper`).
#' @param has_s Logical; `TRUE` for green-myerson / rachlin.
#' @return The (possibly augmented) `tmb_control` list.
#' @keywords internal
.dd_apply_s_bounds <- function(tmb_control, has_s) {
  if (!isTRUE(has_s)) return(tmb_control)
  lo <- tmb_control$lower
  if (is.null(lo) || !("log_s" %in% names(lo))) {
    tmb_control$lower <- c(lo, log_s = .dd_s_log_lower)
  }
  up <- tmb_control$upper
  if (is.null(up) || !("log_s" %in% names(up))) {
    tmb_control$upper <- c(up, log_s = .dd_s_log_upper)
  }
  tmb_control
}


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
#' @param has_s Logical; `TRUE` for the 2-parameter equations.
#' @param n_re Number of random-effect intercepts (`1L` or `2L`). Defaults to
#'   `1L` so existing positional callers stay valid.
#' @param covariance `"pdSymm"` or `"pdDiag"` (2-RE only).
#' @return A list with elements `obj`, `opt`, `nll`, `start_idx`,
#'   `opt_warnings`.
#' @keywords internal
.dd_tmb_multi_start <- function(tmb_data, start_values,
                                 tmb_control, user_specified, verbose,
                                 has_s = FALSE, n_re = 1L,
                                 covariance = "pdSymm") {
  # Impose the phi floor as a log_aux lower bound for sltb (shared helper used by
  # both optimizer paths; tmb_data$family == 0L encodes "sltb").
  family_str <- if (identical(tmb_data$family, 0L)) "sltb" else "gaussian"
  tmb_control <- .dd_apply_phi_floor(tmb_control, family_str)
  tmb_control <- .dd_apply_s_bounds(tmb_control, has_s)
  tmb_control <- .dd_apply_re_bounds(tmb_control, n_re)   # 2-RE degeneracy guard
  random_block <- if (n_re == 2L) "b" else "u"
  map <- .dd_tmb_build_map(has_s, n_re, covariance)
  map <- .dd_tmb_finalize_map(map, start_values, n_re)

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

  # For the 2-parameter equations also perturb the (poorly-identified) log_s
  # axis across the low-k / high-k restarts.
  if (isTRUE(has_s)) {
    start_sets[[2]]$log_s <- log(0.7)
    start_sets[[3]]$log_s <- log(1.4)
  }

  # For a 2-RE fit perturb the (poorly-identified) RE-variance/correlation axes
  # across the restarts, analogous to the log_s perturbation.
  if (n_re == 2L) {
    start_sets[[2]]$log_sd_re <- c(log(0.3), log(0.3))
    start_sets[[3]]$log_sd_re <- c(log(0.8), log(0.8))
    if (identical(covariance, "pdSymm")) {
      start_sets[[2]]$cor_re <- -0.2
      start_sets[[3]]$cor_re <- 0.4
    }
  }

  # Sanity NaN/blowup guard on the FULL fitted log-k predictor (B7) -- no phi
  # rejection (the log_aux lower bound already prevents the degenerate phi->0
  # optimum).
  .is_blowup <- function(opt) .dd_logk_blowup(opt, tmb_data$X)

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
        map        = map,
        random     = random_block,
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
      "All multi-start fits tripped the log-k sanity guard ",
      "(some fitted k -> Inf / non-finite across the design); returning the ",
      "lowest-nll fit. Inspect data for boundary-heavy subjects or degenerate cells.",
      call. = FALSE
    )
  }

  if (verbose >= 1) {
    message(sprintf(
      "  Multi-start: best NLL = %.2f (start set %d of %d)",
      best_result$nll, best_result$start_idx, length(start_sets)
    ))
  }

  # Map-seam assertion: the free-parameter blocks must match (n_re, covariance,
  # has_s) exactly (the aux rename to log_phi/log_sigma_e happens later, so the
  # raw optimizer name log_aux is asserted here).
  expected_free <- .dd_expected_free_params(has_s, n_re, covariance)
  .dd_assert_free_params(names(best_result$opt$par), expected_free,
                         "an IP discounting fit")

  best_result
}


#' Expected free-parameter blocks for a given map configuration
#'
#' The exact set of optimizer parameter blocks left free by the map for a
#' `(has_s, n_re, covariance)` combination. Uses the raw optimizer name
#' `log_aux` (the rename to `log_phi`/`log_sigma_e` happens in
#' [.dd_tmb_extract_estimates()]).
#'
#' @param has_s Logical; `TRUE` for the 2-parameter equations.
#' @param n_re Number of random-effect intercepts (`1L` or `2L`).
#' @param covariance `"pdSymm"` or `"pdDiag"` (2-RE only).
#' @return Character vector of expected free-block names.
#' @keywords internal
.dd_expected_free_params <- function(has_s, n_re, covariance) {
  expected <- c("beta_k", "log_aux")
  if (isTRUE(has_s)) expected <- c(expected, "log_s")
  if (n_re == 2L) {
    expected <- c(expected, "log_sd_re")
    if (identical(covariance, "pdSymm")) expected <- c(expected, "cor_re")
  } else {
    expected <- c(expected, "log_sigma_u")
  }
  expected
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
#' @param has_s Logical; `TRUE` for the 2-parameter equations (green-myerson /
#'   rachlin), where `log_s` is a free coefficient.
#' @param verbose Integer verbosity.
#' @param n_re Number of random-effect intercepts (`1L` or `2L`). Defaults to
#'   `1L` so existing positional callers stay valid.
#' @return list(coefficients, se, sdr, variance_components, u_hat, hessian_pd,
#'   Sigma). `u_hat` is the per-subject random-effect block: a 1-column matrix
#'   of standardized `u` deviates for `n_re == 1`, an `n_subjects x 2` matrix of
#'   standardized `b` deviates for `n_re == 2`. `Sigma` is the fitted 2x2 RE
#'   covariance for `n_re == 2` (NULL otherwise).
#' @keywords internal
.dd_tmb_extract_estimates <- function(obj, opt, n_subjects, family,
                                      has_s = FALSE, verbose = 1, n_re = 1L,
                                      re2_target = 0L) {
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

  # The presence of a FREE log_s in opt$par is authoritative: it must match
  # has_s (a 1-parameter fit maps log_s out; a 2-parameter fit leaves it free).
  free_log_s <- "log_s" %in% par_names
  if (!identical(free_log_s, isTRUE(has_s))) {
    stop("internal: free log_s (", free_log_s, ") does not match has_s (",
         isTRUE(has_s), "); check map threading.", call. = FALSE)
  }

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
    .fill_vector_se("log_s")   # no-op when log_s is mapped/absent
    .fill_vector_se("log_sd_re")  # no-op for a 1-RE fit (mapped out)
    .fill_vector_se("cor_re")     # no-op for a 1-RE / pdDiag fit

    re_summary <- tryCatch(summary(sdr, "random"), error = function(e) NULL)
    re_ncol <- if (n_re == 2L) 2L else 1L
    if (!is.null(re_summary)) {
      re_hat <- matrix(re_summary[, "Estimate"], nrow = n_subjects,
                       ncol = re_ncol)
    } else {
      re_hat <- matrix(0, nrow = n_subjects, ncol = re_ncol)
    }
  } else {
    re_hat <- matrix(0, nrow = n_subjects, ncol = if (n_re == 2L) 2L else 1L)
  }

  variance_components <- NULL
  if (!is.null(sdr)) {
    variance_components <- tryCatch(summary(sdr, "report"), error = function(e) NULL)
  }

  # For a 2-RE fit, rebuild the fitted RE covariance from the fitted
  # COEFFICIENTS (log_sd_re / cor_re), not the sdreport "report" block. The
  # coefficients are always available on a converged fit, so Sigma is non-NULL
  # even if sdreport fails -- matching the kernel's Sigma = diag(sd) R diag(sd).
  Sigma <- NULL
  if (n_re == 2L) {
    sd_re_hat <- exp(unname(coefficients[names(coefficients) == "log_sd_re"]))
    cor_hat   <- if ("cor_re" %in% names(coefficients)) {
      unname(coefficients[["cor_re"]])
    } else {
      0
    }
    rho_hat <- tanh(cor_hat)            # pdDiag: cor_re mapped out -> rho = 0
    re2_name <- if (identical(as.integer(re2_target), 1L)) "s" else "phi"
    Sigma <- matrix(c(sd_re_hat[1]^2, rho_hat * sd_re_hat[1] * sd_re_hat[2],
                      rho_hat * sd_re_hat[1] * sd_re_hat[2], sd_re_hat[2]^2), 2L,
                    dimnames = list(c("k", re2_name), c("k", re2_name)))
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
    u_hat = re_hat,
    Sigma = Sigma,
    hessian_pd = hessian_pd
  )
}


#' Compute subject-specific discounting parameters
#'
#' Reconstructs each subject's `k_i = exp(eta_i)` where
#' `eta_i = X_i %*% beta_k + sigma_u * u_i`, exactly matching the C++ template's
#' non-centered predictor (`log_k_i = X.row(i)*beta_k + sigma_u * u(subj,0)`).
#' `X_i` is subject `i`'s fixed-effect design ROW. For between-subject
#' factors/covariates the design is constant within a subject, so the first
#' design-matrix row for each subject is used. This makes per-subject `k`
#' factor/covariate-correct: subjects in non-reference groups pick up their
#' group's `beta_k` contribution rather than the reference-group intercept.
#'
#' For a 1-RE fit the auxiliary scalar `phi` is population-level, so it is not a
#' subject-level parameter (the returned frame is `id, u_i, k`). For a φ-target
#' 2-RE fit (`k + phi ~ 1`) each subject's natural-scale `re = L * b_hat` is
#' reconstructed (`L = chol(Sigma)` lower-triangular; the `sdreport` "random"
#' block holds the STANDARDIZED `b_hat`, not `re`), giving `id, re_k, re_phi,
#' k, phi` with a per-subject `phi_i` floor matching the kernel clamp. For an
#' s-target 2-RE fit (`k + s ~ 1`) the same reconstruction gives `id, re_k,
#' re_s, k, s` with `s_i = exp(log_s + re_s)` clamped to `[0.05, 20]`.
#'
#' @param coefficients Named coefficient vector (with `beta_k`, `log_sigma_u`
#'   for 1-RE, and `log_phi` or `log_sigma_e`).
#' @param u_hat Matrix of per-subject random effects: `n_subjects x 1`
#'   standardized `u` deviates (1-RE) or `n_subjects x 2` standardized `b`
#'   deviates (2-RE).
#' @param subject_levels Character vector of subject ids (length n_subjects).
#' @param design_X Fixed-effect design matrix (rows aligned with the prepared
#'   arrays); columns correspond 1:1 with the `beta_k` entries.
#' @param subject_id Integer 0-indexed subject map (length `n_obs`), aligned
#'   row-for-row with `design_X`. Used to locate each subject's first design row.
#' @param equation One of "mazur", "exponential" (reserved; k is equation-free).
#' @param family One of "sltb", "gaussian".
#' @param n_re Number of random-effect intercepts (`1L` or `2L`). Defaults to
#'   `1L` so existing positional callers stay valid.
#' @param log_aux_name Name of the auxiliary (precision) coefficient in
#'   `coefficients` (`"log_phi"` for sltb). Used only for `n_re == 2`.
#' @param Sigma Fitted 2x2 RE covariance (from [.dd_tmb_extract_estimates()]);
#'   required for `n_re == 2`.
#' @param re2_target Integer flag: `0L` for a φ-target 2-RE fit, `1L` for an
#'   s-target fit. Ignored when `n_re == 1L`. Defaults to `0L`.
#' @return data.frame: `id, u_i, k` (1-RE), `id, re_k, re_phi, k, phi`
#'   (φ-target 2-RE), or `id, re_k, re_s, k, s` (s-target 2-RE).
#' @note Subject-level `k` assumes between-subject predictors: the first design
#'   row per subject defines that subject's fixed-effect contribution. A
#'   within-subject-varying covariate would make a single per-subject `k`
#'   ill-defined; in that case the first observed row is used.
#' @keywords internal
.dd_tmb_compute_subject_pars <- function(coefficients, u_hat, subject_levels,
                                         design_X, subject_id,
                                         equation, family, n_re = 1L,
                                         log_aux_name = "log_phi",
                                         Sigma = NULL, re2_target = 0L) {
  beta_k <- unname(coefficients[names(coefficients) == "beta_k"])
  n_subjects <- length(subject_levels)

  design_X <- as.matrix(design_X)

  # First design row per subject (subject_id is 0-indexed: subject s -> s - 1).
  # Between-subject predictors are constant within a subject, so the first row
  # carries that subject's fixed-effect contribution; for any within-subject
  # variation the first observed row is used (documented in @note).
  xbeta_i <- numeric(n_subjects)
  for (s in seq_len(n_subjects)) {
    first_row <- which(subject_id == (s - 1L))[1L]
    xbeta_i[s] <- sum(design_X[first_row, ] * beta_k)
  }

  if (n_re == 2L) {
    if (is.null(Sigma)) stop("internal: 2-RE Sigma unavailable", call. = FALSE)
    b_hat <- u_hat                       # n_subjects x 2 standardized deviates
    L  <- t(chol(Sigma))                 # lower-triangular Cholesky of fitted Sigma
    re <- t(L %*% t(b_hat))             # natural-scale (re_k, re_2)
    re_k <- as.numeric(re[, 1L])
    log_k_i <- xbeta_i + re_k
    if (identical(as.integer(re2_target), 1L)) {
      # s-target: second RE is on log s; phi/sigma_e stays population.
      re_s  <- as.numeric(re[, 2L])
      log_s <- unname(coefficients[["log_s"]])
      s_i   <- pmin(pmax(exp(log_s + re_s), 0.05), 20)   # same clamp as the kernel
      data.frame(id = subject_levels, re_k = re_k, re_s = re_s,
                 k = exp(log_k_i), s = s_i, stringsAsFactors = FALSE)
    } else {
      re_phi  <- as.numeric(re[, 2L])
      log_aux <- unname(coefficients[[log_aux_name]])
      phi_i   <- pmax(exp(log_aux + re_phi), 0.1)        # per-subject phi floor
      data.frame(id = subject_levels, re_k = re_k, re_phi = re_phi,
                 k = exp(log_k_i), phi = phi_i, stringsAsFactors = FALSE)
    }
  } else {
    sigma_u <- exp(unname(coefficients[["log_sigma_u"]]))
    u_i <- as.numeric(u_hat[, 1L])
    log_k_i <- xbeta_i + sigma_u * u_i
    data.frame(id = subject_levels, u_i = u_i, k = exp(log_k_i),
               stringsAsFactors = FALSE)
  }
}


#' Warn when a declared predictor is not constant within a subject (B4)
#'
#' `fit_dd_tmb()` treats factors/covariates as BETWEEN-subject, and
#' `.dd_tmb_compute_subject_pars()` reconstructs each subject's `k` from that
#' subject's FIRST design row. A predictor that varies within a subject still
#' yields a valid row-varying fit, but the per-subject `k` / `ranef()` summary is
#' then approximate, so this warns (it does not abort) naming the offending
#' column(s). Numeric columns use a tolerance so floating noise does not warn.
#'
#' @param data The single complete-cased model frame (canonical `id`/`x`/`y` plus
#'   any retained factor/covariate columns).
#' @param extra_cols Character vector of declared factor/covariate column names.
#' @param id_col Subject id column (canonical `"id"`).
#' @param tol Numeric tolerance for "constant within id" on numeric columns.
#' @return Invisibly, the character vector of within-subject-varying columns.
#' @keywords internal
#' @noRd
.dd_check_between_subject <- function(data, extra_cols, id_col = "id",
                                      tol = 1e-8) {
  varying <- character(0)
  affected <- character(0)
  for (col in intersect(unique(extra_cols), names(data))) {
    not_const <- tapply(data[[col]], data[[id_col]], function(v) {
      v <- v[!is.na(v)]
      if (length(v) <= 1L) return(FALSE)
      if (is.numeric(v)) (max(v) - min(v)) > tol else length(unique(v)) > 1L
    })
    if (any(not_const, na.rm = TRUE)) {
      varying <- c(varying, col)
      affected <- union(affected, names(not_const)[which(not_const)])
    }
  }
  if (length(varying) > 0L) {
    n_affected <- length(affected)
    cli::cli_warn(c(
      "!" = "{cli::qty(varying)}Predictor{?s} {.val {varying}} {?is/are} not \\
             constant within {.field {id_col}} ({cli::qty(n_affected)}{n_affected} \\
             subject{?s} affected).",
      "i" = "{.fn fit_dd_tmb} treats predictors as between-subject; subject-level \\
             {.code k} / {.fn ranef} use each subject's FIRST design row, so they \\
             are approximate for within-subject-varying predictors."
    ))
  }
  invisible(varying)
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
#' @param equation One of `"mazur"`, `"exponential"`, `"green-myerson"`, or
#'   `"rachlin"`. The two 2-parameter (hyperboloid) forms add a single
#'   population nonlinearity exponent `s` (estimated on the log scale) and
#'   reduce to `"mazur"` at `s = 1`.
#' @param family Observation family: `"sltb"` (default) or `"gaussian"`.
#' @param random_effects RE formula: `k ~ 1` (single random intercept on
#'   `log k`), `k + phi ~ 1` (a joint 2-D random intercept on
#'   `(log k, log phi)`, SLT-beta only), or `k + s ~ 1` (a joint 2-D random
#'   intercept on `(log k, log s)`, Green-Myerson and Rachlin only).
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
#' @param covariance_structure Covariance for a 2-D random effect
#'   (`k + phi ~ 1` or `k + s ~ 1`): `"pdSymm"` (default; correlated random
#'   intercepts) or `"pdDiag"` (independent, correlation fixed at 0). Ignored
#'   for `k ~ 1`.
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
#'     \item{subject_pars}{Data frame of subject-level parameters. For a 1-RE fit
#'       (`k ~ 1`) the columns are `id, u_i, k`; for a φ-target 2-RE fit
#'       (`k + phi ~ 1`) they are `id, re_k, re_phi, k, phi`; for an s-target
#'       2-RE fit (`k + s ~ 1`, GM/Rachlin) they are `id, re_k, re_s, k, s`
#'       where `s` is clamped to `[0.05, 20]`.}
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
                       equation = c("mazur", "exponential",
                                    "green-myerson", "rachlin"),
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
                       covariance_structure = c("pdSymm", "pdDiag"),
                       ...) {
  cl <- match.call()
  equation <- match.arg(equation)
  family <- match.arg(family)
  covariance_structure <- match.arg(covariance_structure)
  response_scale <- match.arg(response_scale)
  has_s <- equation %in% c("green-myerson", "rachlin")

  # The union of factor + covariate columns is the set of EXTRA modeling
  # columns that must travel with id/x/y through one complete-case pass
  # (row-coherence). The validator retains them, prepare_data complete-cases
  # them ONCE, and the design + TMB arrays are derived from that single frame.
  extra_cols <- unique(c(factors, continuous_covariates))

  # B5/B6: declared predictors must not BE the id/delay/response columns -- e.g.
  # factors = "x" would otherwise coerce the delay to a factor below.
  # (.dd_validate_ip separately reserves the canonical id/x/y output names for
  # the remapped-role case.)
  role_collision <- intersect(extra_cols, c(id_var, x_var, y_var))
  if (length(role_collision) > 0L) {
    cli::cli_abort(c(
      "{.arg factors}/{.arg continuous_covariates} cannot be the id, delay, or \\
       response column{?s}: {.val {role_collision}}.",
      "i" = "These name the subject/delay/response, not between-subject predictors."
    ))
  }

  # 1. Validate + coerce/clamp (warns; errors on missing/NA y). Retains the
  #    factor/covariate columns alongside canonical id/x/y.
  validated <- .dd_validate_ip(data, y_var = y_var, x_var = x_var,
                               id_var = id_var, ll = ll,
                               extra_cols = extra_cols,
                               response_scale = response_scale)
  long <- validated$data            # canonical id/x/y + retained extras
  coercion_info <- validated$coercion_info

  # 2. Random-effects normalization: `k ~ 1` (1-RE), or `k + phi ~ 1` / `k + s ~ 1` (2-RE).
  re_norm <- .dd_normalize_re(random_effects, covariance_structure, data = long)
  n_random_effects <- re_norm$blocks[[1]]$dim
  re_has_phi <- "phi" %in% re_norm$blocks[[1]]$param
  if (re_has_phi && family != "sltb") {
    cli::cli_abort(c(
      "A subject-random {.field phi} ({.code k + phi ~ 1}) requires \\
       {.code family = \"sltb\"}.",
      "i" = "The Gaussian family has no precision parameter to put a random \\
             effect on."
    ))
  }
  re_has_s <- "s" %in% re_norm$blocks[[1]]$param
  if (re_has_s && !has_s) {
    cli::cli_abort(c(
      "A subject-random {.field s} ({.code k + s ~ 1}) requires \\
       {.code equation = \"green-myerson\"} or {.code \"rachlin\"}.",
      "i" = "The Mazur and exponential equations have no {.field s} curvature \\
             parameter to put a random effect on."
    ))
  }
  re2_target <- if (re_has_s) 1L else 0L
  re_cov <- re_norm$blocks[[1]]$pdmat_class      # "pdDiag" or "pdSymm"

  # 3. Prepare data: ONE complete-case pass over id/x/y + extra_cols, building
  #    the 0-indexed subject_id. prepared$data is the single filtered model
  #    frame (id/x/y + retained extras) that everything else derives from.
  prepared <- .dd_tmb_prepare_data(long, y_var = "y", x_var = "x",
                                   id_var = "id", extra_cols = extra_cols)

  # 3b. B5: coerce declared factors to factor (a numeric column passed as a
  #     factor would otherwise fit as a single continuous slope, not level
  #     dummies), and require continuous covariates to be numeric AND finite
  #     (Inf survives the complete-case filter). Operate on prepared$data so the
  #     design and the stored fit$data share the coerced types.
  for (f in factors) {
    if (f %in% names(prepared$data) && !is.factor(prepared$data[[f]])) {
      prepared$data[[f]] <- as.factor(prepared$data[[f]])
    }
  }
  for (cv in continuous_covariates) {
    if (cv %in% names(prepared$data)) {
      v <- prepared$data[[cv]]
      if (!is.numeric(v) || any(!is.finite(v))) {
        cli::cli_abort(c(
          "Continuous covariate {.val {cv}} must be numeric and finite.",
          "i" = "Found {.val {if (!is.numeric(v)) class(v)[1] else 'non-finite value(s)'}}."
        ))
      }
    }
  }

  # 3c. B4: between-subject contract. A predictor that varies WITHIN a subject
  #     still fits (row-varying GLM) but makes subject_pars$k / ranef() reflect
  #     only the first design row, so warn (not abort) naming the column(s).
  .dd_check_between_subject(prepared$data, extra_cols)

  # 4. Fixed-effect design for log k, built on the SAME filtered frame
  #    (prepared$data) so X's rows align 1:1 with the prepared y/x/subject_id.
  design <- .dd_tmb_build_design(prepared$data, factors = factors,
                                 factor_interaction = factor_interaction,
                                 continuous_covariates = continuous_covariates)

  # Record only the factors actually placed in the design RHS (build_fixed_rhs
  # drops single-level factors). param_info$factors must reflect the design.
  placed_factors <- intersect(factors, all.vars(design$rhs))
  if (length(placed_factors) == 0L) placed_factors <- NULL

  # 5. TMB data + default starts.
  tmb_data <- .dd_tmb_build_tmb_data(prepared, design, equation, family,
                                     n_re = n_random_effects,
                                     re2_target = re2_target)
  default_starts <- .dd_tmb_default_starts(prepared, design, family, equation)
  if (!is.null(start_values)) {
    # Public alias: users pass `s` (natural); convert to the optimizer's log_s.
    if (!is.null(start_values$s)) {
      if (!is.null(start_values$log_s)) {
        stop("start_values: supply either 's' or 'log_s', not both.", call. = FALSE)
      }
      start_values$log_s <- log(start_values$s)
      start_values$s <- NULL
    }
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
                                  user_specified, verbose, has_s = has_s,
                                  n_re = n_random_effects, covariance = re_cov)
    obj <- result$obj
    opt <- result$opt
    opt_warnings <- result$opt_warnings %||% character(0)
  } else {
    # Single-path fit still needs the SLT-beta phi floor (B3): apply the same
    # shared log_aux lower bound the multi-start path uses, so multi_start=FALSE
    # cannot collapse into the degenerate phi->0 / k->inf optimum.
    tmb_control <- .dd_apply_phi_floor(tmb_control, family)
    tmb_control <- .dd_apply_s_bounds(tmb_control, has_s)
    tmb_control <- .dd_apply_re_bounds(tmb_control, n_random_effects)
    random_block <- if (n_random_effects == 2L) "b" else "u"
    single_map <- .dd_tmb_build_map(has_s, n_random_effects, re_cov)
    single_map <- .dd_tmb_finalize_map(single_map, default_starts,
                                       n_random_effects)
    obj <- TMB::MakeADFun(tmb_data, default_starts, random = random_block,
                          DLL = "beezdiscounting", silent = verbose < 2,
                          map = single_map)
    opt_res <- .dd_tmb_run_optimizer(obj, obj$par, tmb_control,
                                     user_specified, verbose)
    opt <- opt_res$opt
    opt_warnings <- opt_res$warnings
    .dd_assert_free_params(
      names(opt$par),
      .dd_expected_free_params(has_s, n_random_effects, re_cov),
      "an IP discounting fit")
    # B7: the single-start path has no multi-start selection guard, so warn if the
    # fitted log-k predictor blew up (k -> Inf / non-finite).
    if (.dd_logk_blowup(opt, tmb_data$X)) {
      warning(
        "The fitted log-k linear predictor exceeds the sanity bound ",
        "(k -> Inf / non-finite); inspect data for boundary-heavy subjects or ",
        "degenerate cells.",
        call. = FALSE
      )
    }
  }

  converged <- isTRUE(opt$convergence == 0)
  try(obj$fn(opt$par), silent = TRUE)

  # 8. Extract estimates (sdreport, pdHess gate, log_aux rename).
  estimates <- .dd_tmb_extract_estimates(obj, opt,
                                         n_subjects = prepared$n_subjects,
                                         family = family, has_s = has_s,
                                         verbose = verbose,
                                         n_re = n_random_effects,
                                         re2_target = re2_target)

  # The aux scalar's renamed coefficient name (mirrors the rename in
  # .dd_tmb_extract_estimates) -- needed for the 2-RE per-subject phi.
  aux_name <- if (identical(family, "gaussian")) "log_sigma_e" else "log_phi"

  # 9. Subject-specific parameters. 1-RE: id/u_i/k (phi is population-level).
  #    2-RE: id/re_k/re_phi/k/phi (per-subject phi). Factor/covariate-correct:
  #    each subject's k uses that subject's design row (X_i %*% beta_k + re_k),
  #    not the reference-group intercept (B1).
  subject_pars <- .dd_tmb_compute_subject_pars(
    coefficients = estimates$coefficients,
    u_hat = estimates$u_hat,
    subject_levels = prepared$subject_levels,
    design_X = design$X,
    subject_id = prepared$subject_id,
    equation = equation,
    family = family,
    n_re = n_random_effects,
    log_aux_name = aux_name,
    Sigma = estimates$Sigma,
    re2_target = re2_target
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
      # Fitted 2x2 RE covariance (k, phi) for a 2-RE fit; NULL for 1-RE. Read by
      # VarCorr.beezdiscounting_tmb and the per-subject re reconstruction.
      Sigma = estimates$Sigma,
      hessian_pd = estimates$hessian_pd,
      param_info = list(
        equation = equation,
        has_s = has_s,
        family = family,
        has_phi = has_phi,
        n_obs = prepared$n_obs,
        n_subjects = prepared$n_subjects,
        n_random_effects = n_random_effects,
        covariance_structure = re_cov,
        re2_target = re2_target,
        re_aux_param = if (re_has_s) "s" else if (re_has_phi) "phi" else NA_character_,
        subject_levels = prepared$subject_levels,
        # B1: store the CANONICAL column names. The validator renames the data to
        # id/x/y, and fit$data (the default newdata for predict/fitted/residuals/
        # augment) is always canonical, so the post-fit methods must index it by
        # the canonical names -- not the caller's original y_var/x_var/id_var,
        # which would not exist on fit$data and broke every response-path method.
        # The user's original names are kept under user_vars for reference/display.
        id_var = "id",
        x_var = "x",
        y_var = "y",
        user_vars = list(id = id_var, x = x_var, y = y_var),
        # Factor metadata persisted so predict()/emmeans can rebuild a
        # column-aligned design via build_fixed_rhs + stored contrasts (B3).
        # Only factors actually placed in the design are recorded (single-level
        # factors are dropped by build_fixed_rhs).
        factors = placed_factors,
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
      # SE-availability requires BOTH a successful sdreport AND a positive-
      # definite Hessian (R1): a non-PD Hessian yields untrustworthy SEs, so
      # SE-consuming methods must not present them as reliable.
      se_available = !is.null(estimates$sdr) && isTRUE(estimates$hessian_pd),
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
