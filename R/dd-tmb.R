#' @useDynLib beezdiscounting, .registration = TRUE
NULL

# ==============================================================================
# Internal Helper Functions for TMB Mixed-Effects Discounting Models
# ==============================================================================

#' Prepare data for the TMB mixed-effects discounting model
#'
#' Operates on the FULL model frame from [.dd_validate_ip()] (canonical
#' `id`/`x`/`y` PLUS any retained factor/covariate columns). It complete-cases
#' that frame ONCE over all modeling columns (`id`, `x`, `y`, and `extra_cols`),
#' then derives the parallel arrays `y`/`x`/`subject_id` from the same filtered
#' rows so the later design matrix, TMB arrays, and `fit$data` all share one
#' row order (row-coherence; mirrors beezdemand's `data_for_design`). Builds a
#' 0-indexed `subject_id` aligned to `subject_levels` (the C++ template indexes
#' `u(subject_id, 0)` from 0).
#'
#' @param data Data frame already coerced/clamped by [.dd_validate_ip()],
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
#'   `contrasts` are stored so [.dd_build_emm_ref_grid()] can rebuild a
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
