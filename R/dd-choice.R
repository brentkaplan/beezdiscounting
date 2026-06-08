# ==============================================================================
# Choice family (Family 2) — trial-level SS-vs-LL binomial GLMM (structural).
# ==============================================================================

#' Validate + coerce a trial-level SS-vs-LL choice data frame
#'
#' Remaps caller columns to the canonical `id, ss_amount, ll_amount, delay,
#' choice` frame, coerces types, and loudly rejects malformed data: `choice`
#' must be binary `{0,1}`; `ss_amount`/`ll_amount` finite & > 0; `delay` finite
#' & > 0 (SS is immediate, delay 0, this slice). `ll_amount <= ss_amount` warns
#' (not an error — some designs include catch trials). Declared `extra_cols`
#' (between-subject factors/covariates for the log-k design) are retained
#' verbatim (class/levels preserved).
#'
#' @param data Data frame with the choice columns.
#' @param id_var,ss_var,ll_var,delay_var,choice_var Column names.
#' @param extra_cols Optional character vector of additional columns to retain
#'   (between-subject factors/covariates).
#' @return list(data = canonical frame (+ extras), coercion_info = list()).
#' @keywords internal
#' @noRd
.dd_validate_choice <- function(data, id_var = "id", ss_var = "ss_amount",
                                ll_var = "ll_amount", delay_var = "delay",
                                choice_var = "choice", extra_cols = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  needed <- c(id_var, ss_var, ll_var, delay_var, choice_var)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(shQuote(missing_cols), collapse = ", "),
         ". Set the *_var arguments to match your data.", call. = FALSE)
  }
  long <- data.frame(
    id = as.character(data[[id_var]]),
    ss_amount = suppressWarnings(as.numeric(data[[ss_var]])),
    ll_amount = suppressWarnings(as.numeric(data[[ll_var]])),
    delay = suppressWarnings(as.numeric(data[[delay_var]])),
    choice = suppressWarnings(as.numeric(data[[choice_var]])),
    stringsAsFactors = FALSE
  )

  # choice strictly binary {0,1} (NA tolerated here; dropped in prepare).
  ch <- long$choice
  bad_choice <- !is.na(ch) & !(ch %in% c(0, 1))
  if (any(bad_choice)) {
    stop(sprintf("`%s` must be binary 0/1 (1 = chose LL); %d value(s) are not.",
                 choice_var, sum(bad_choice)), call. = FALSE)
  }

  # amounts finite & > 0.
  for (pair in list(c("ss_amount", ss_var), c("ll_amount", ll_var))) {
    v <- long[[pair[[1]]]]
    if (any(!is.na(v) & (!is.finite(v) | v <= 0))) {
      stop(sprintf("`%s` must be finite and positive.", pair[[2]]), call. = FALSE)
    }
  }
  # delay finite & > 0.
  d <- long$delay
  if (any(!is.na(d) & !is.finite(d))) {
    stop(sprintf("`%s` must be finite.", delay_var), call. = FALSE)
  }
  if (any(!is.na(d) & d <= 0)) {
    stop(sprintf("`%s` must be > 0 (SS is immediate; LL delay positive).",
                 delay_var), call. = FALSE)
  }
  # ll > ss expected (warn).
  bad_order <- !is.na(long$ll_amount) & !is.na(long$ss_amount) &
    long$ll_amount <= long$ss_amount
  if (any(bad_order)) {
    warning(sprintf(paste0("%d row(s) have ll_amount <= ss_amount; expected the ",
                           "larger-later reward to exceed the smaller-sooner one."),
                    sum(bad_order)), call. = FALSE)
  }

  # Retain declared extra columns (between-subject design) verbatim.
  extra_cols <- setdiff(unique(extra_cols),
                        c(id_var, ss_var, ll_var, delay_var, choice_var))
  if (length(extra_cols)) {
    missing_extra <- setdiff(extra_cols, names(data))
    if (length(missing_extra)) {
      stop("Factor/covariate column(s) not found: ",
           paste(shQuote(missing_extra), collapse = ", "), ".", call. = FALSE)
    }
    collide <- intersect(extra_cols,
                         c("id", "ss_amount", "ll_amount", "delay", "choice"))
    if (length(collide)) {
      stop("Factor/covariate column(s) collide with reserved canonical names: ",
           paste(shQuote(collide), collapse = ", "), ".", call. = FALSE)
    }
    for (col in extra_cols) long[[col]] <- data[[col]]  # preserve class/levels
  }

  list(data = long, coercion_info = list(n_rows_in = nrow(data)))
}


#' Prepare validated choice data for TMB (0-indexed subjects, row-coherence)
#'
#' Complete-cases over all modeling columns (canonical + declared extras) ONCE,
#' builds a 0-indexed `subject_id` aligned to sorted `subject_levels` (index by
#' name, never positionally), and returns parallel arrays + the single filtered
#' frame. Retained factor extras are `droplevels()`-ed after row filtering.
#'
#' @param data Canonical frame (+ extras) from [.dd_validate_choice()].
#' @param extra_cols Optional character vector of extra columns to keep.
#' @return list(ss_amount, ll_amount, delay, choice, subject_id, subject_levels,
#'   n_subjects, n_obs, data).
#' @keywords internal
#' @noRd
.dd_choice_prepare_data <- function(data, extra_cols = NULL) {
  canonical <- c("id", "ss_amount", "ll_amount", "delay", "choice")
  model_cols <- intersect(unique(c(canonical, extra_cols)), names(data))
  keep <- stats::complete.cases(data[, model_cols, drop = FALSE])
  frame <- data[keep, , drop = FALSE]
  if (nrow(frame) == 0L) {
    stop("No complete cases remain after dropping NA rows.", call. = FALSE)
  }
  subject_levels <- sort(unique(frame$id))
  subject_map <- stats::setNames(seq_along(subject_levels) - 1L, subject_levels)
  subject_id <- as.integer(subject_map[frame$id])
  extras <- setdiff(intersect(extra_cols, names(frame)), canonical)
  cleaned <- frame[, canonical, drop = FALSE]
  for (col in extras) {
    v <- frame[[col]]
    if (is.factor(v)) v <- droplevels(v)
    cleaned[[col]] <- v
  }
  list(
    ss_amount = cleaned$ss_amount, ll_amount = cleaned$ll_amount,
    delay = cleaned$delay, choice = cleaned$choice,
    subject_id = subject_id, subject_levels = subject_levels,
    n_subjects = length(subject_levels), n_obs = nrow(cleaned),
    data = cleaned
  )
}


#' Discount function value for the structural choice model (mazur/exponential)
#'
#' This slice supports the two 1-parameter equations only (GM/Rachlin deferred).
#' SS is immediate, so `D(k, 0) = 1` and `V_SS = ss_amount`.
#' @keywords internal
#' @noRd
.dd_choice_D <- function(k, delay, equation) {
  switch(equation,
    mazur = 1 / (1 + k * delay),
    exponential = exp(-k * delay),
    stop("unsupported equation '", equation,
         "' (structural choice: mazur or exponential)", call. = FALSE)
  )
}

#' Structural choice linear predictor (logit scale), scale-invariant form
#'
#' `eta = beta0 + gamma * ((ll/ss) * D(k, delay) - 1)`, i.e.
#' `beta0 + gamma * (V_LL - V_SS) / V_SS` with `V_SS = ss_amount`,
#' `V_LL = ll_amount * D(k, delay)`. Magnitude enters only through the `ll/ss`
#' ratio (scale-invariant). With `beta0 = 0`, `eta = 0` <=> `V_LL = V_SS`
#' (the indifference point), so `k` is the classical discount rate.
#'
#' @param k Numeric vector of per-row discount rates (or scalar recycled).
#' @param ss_amount,ll_amount,delay Numeric per-row trial values.
#' @param equation `"mazur"` or `"exponential"`.
#' @param gamma Choice sensitivity (> 0).
#' @param beta0 Choice-bias intercept (0 when `intercept = FALSE`).
#' @return Numeric vector of logit-scale linear predictors.
#' @keywords internal
#' @noRd
.dd_choice_structural_eta <- function(k, ss_amount, ll_amount, delay, equation,
                                      gamma, beta0 = 0) {
  D <- .dd_choice_D(k, delay, equation)
  beta0 + gamma * ((ll_amount / ss_amount) * D - 1)
}


# ==============================================================================
# Structural choice fit pipeline (direct-k binomial GLMM via TMB).
# ==============================================================================

#' TMB `map` for the structural choice model: fix beta0 when intercept is off
#' @keywords internal
#' @noRd
.dd_choice_build_map <- function(intercept) {
  if (isTRUE(intercept)) return(NULL)
  list(beta0 = factor(NA))
}

#' Default starts for the structural choice model
#' @keywords internal
#' @noRd
.dd_choice_default_starts <- function(prepared, design) {
  p <- ncol(design$X)
  beta_k <- rep(0, p)
  beta_k[1] <- log(0.01)                  # generic mid-range k start
  list(
    beta_k = beta_k, log_sigma_u = log(0.5),
    log_gamma = log(1), beta0 = 0,
    u = matrix(0, nrow = prepared$n_subjects, ncol = 1L)
  )
}

#' Extract estimates from a structural choice fit (sdreport + SEs)
#' @keywords internal
#' @noRd
.dd_choice_extract_estimates <- function(obj, opt, n_subjects, intercept,
                                         verbose = 1) {
  sdr <- tryCatch(TMB::sdreport(obj), error = function(e) NULL)
  hessian_pd <- if (!is.null(sdr)) isTRUE(sdr$pdHess) else NA
  if (isTRUE(verbose >= 1)) {
    if (is.null(sdr)) {
      cli::cli_warn(c("Standard errors unavailable: {.fn TMB::sdreport} failed.",
                      "i" = "Fixed-effect SEs/CIs will be {.val NA}."))
    } else if (!isTRUE(hessian_pd)) {
      cli::cli_warn(c("Standard errors may be unreliable: the Hessian is not \\
                       positive-definite.",
                      "i" = "Fixed-effect SEs/CIs will be {.val NA}."))
    }
  }
  par_full <- opt$par
  par_names <- names(par_full)
  free_beta0 <- "beta0" %in% par_names
  if (!identical(free_beta0, isTRUE(intercept))) {
    stop("internal: free beta0 (", free_beta0, ") != intercept (",
         isTRUE(intercept), "); check map threading.", call. = FALSE)
  }
  coefficients <- par_full
  se_vec <- stats::setNames(rep(NA_real_, length(par_full)), par_names)
  if (!is.null(sdr)) {
    fixed_summary <- summary(sdr, "fixed")
    fill <- function(name) {
      idx <- which(par_names == name)
      rows <- fixed_summary[rownames(fixed_summary) == name, , drop = FALSE]
      if (length(idx) && nrow(rows) == length(idx)) se_vec[idx] <<- rows[, "Std. Error"]
    }
    fill("beta_k")
    fill("log_sigma_u")
    fill("log_gamma")
    fill("beta0")
    re_summary <- tryCatch(summary(sdr, "random"), error = function(e) NULL)
    u_hat <- if (!is.null(re_summary)) {
      matrix(re_summary[, "Estimate"], nrow = n_subjects, ncol = 1L)
    } else {
      matrix(0, n_subjects, 1L)
    }
  } else {
    u_hat <- matrix(0, n_subjects, 1L)
  }
  list(coefficients = coefficients, se = se_vec, sdr = sdr,
       u_hat = u_hat, hessian_pd = hessian_pd)
}

#' Fit a structural SS-vs-LL choice model (binomial GLMM) via TMB
#'
#' Estimates the discount rate `k` directly from trial-level binary choices via
#' the scale-invariant relative value comparison
#' `logit P(LL) = beta0 + gamma * ((ll/ss) * D(k, delay) - 1)`,
#' `k = exp(X beta_k + sigma_u u)`. Shares the IP family's `k`/emmeans contract.
#'
#' @param data Trial-level data frame (see the `*_var` args).
#' @param mode `"structural"` (this release). `"descriptive"` errors (Plan B).
#' @param id_var,ss_var,ll_var,delay_var,choice_var Column names.
#' @param equation `"mazur"` or `"exponential"`.
#' @param intercept Logical; include the choice-bias `beta0` (default `FALSE`).
#' @param factors,factor_interaction,continuous_covariates Between-subject design
#'   on `log k` (same semantics as [fit_dd_tmb()]).
#' @param start_values,tmb_control,multi_start,verbose,... As in [fit_dd_tmb()].
#' @return An object of class `beezdiscounting_choice`.
#' @export
fit_dd_choice <- function(data, mode = c("structural", "descriptive"),
                          id_var = "id", ss_var = "ss_amount",
                          ll_var = "ll_amount", delay_var = "delay",
                          choice_var = "choice",
                          equation = c("mazur", "exponential"),
                          intercept = FALSE,
                          factors = NULL, factor_interaction = FALSE,
                          continuous_covariates = NULL,
                          start_values = NULL,
                          tmb_control = list(iter_max = 1000, eval_max = 2000),
                          multi_start = TRUE, verbose = 1, ...) {
  cl <- match.call()
  mode <- match.arg(mode)
  equation <- match.arg(equation)
  if (mode == "descriptive") {
    cli::cli_abort(c(
      "{.code mode = \"descriptive\"} is not yet implemented.",
      "i" = "The descriptive (Young 2018) model + random-slope machinery ship \\
             in a later release (Plan B)."
    ))
  }

  # R1: validate + prepare retaining factor/covariate columns for the log-k design
  extra_cols <- unique(c(factors, continuous_covariates))
  validated <- .dd_validate_choice(data, id_var = id_var, ss_var = ss_var,
                                   ll_var = ll_var, delay_var = delay_var,
                                   choice_var = choice_var, extra_cols = extra_cols)
  prepared <- .dd_choice_prepare_data(validated$data, extra_cols = extra_cols)
  for (f in factors) {
    if (f %in% names(prepared$data) && !is.factor(prepared$data[[f]])) {
      prepared$data[[f]] <- as.factor(prepared$data[[f]])
    }
  }
  for (cv in continuous_covariates) {
    if (cv %in% names(prepared$data)) {
      v <- prepared$data[[cv]]
      if (!is.numeric(v) || any(!is.finite(v))) {
        cli::cli_abort("Continuous covariate {.val {cv}} must be numeric and finite.")
      }
    }
  }
  .dd_check_between_subject(prepared$data, extra_cols)
  design <- .dd_tmb_build_design(prepared$data, factors = factors,
                                 factor_interaction = factor_interaction,
                                 continuous_covariates = continuous_covariates)

  tmb_data <- list(
    model = "ChoiceDiscounting", mode = 0L,
    eqn_type = if (equation == "mazur") 0L else 1L,
    has_intercept = as.integer(isTRUE(intercept)),
    choice = as.numeric(prepared$choice),
    ss_amount = as.numeric(prepared$ss_amount),
    ll_amount = as.numeric(prepared$ll_amount),
    delay = as.numeric(prepared$delay),
    subject_id = as.integer(prepared$subject_id),
    X = as.matrix(design$X),
    n_obs = as.integer(prepared$n_obs),
    n_subjects = as.integer(prepared$n_subjects)
  )
  starts <- .dd_choice_default_starts(prepared, design)
  if (!is.null(start_values)) {
    for (nm in names(start_values)) if (nm %in% names(starts)) starts[[nm]] <- start_values[[nm]]
  }
  map <- .dd_choice_build_map(intercept)

  default_control <- list(iter_max = 1000, eval_max = 2000, optimizer = "nlminb",
                          rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
  user_specified <- names(tmb_control)
  tmb_control <- utils::modifyList(default_control, tmb_control)

  perturb_beta_k <- function(delta) {
    b <- starts$beta_k
    b[1] <- b[1] + delta
    b
  }
  start_sets <- list(
    starts,
    utils::modifyList(starts, list(beta_k = perturb_beta_k(-1.5),
                                   log_sigma_u = log(0.3))),
    utils::modifyList(starts, list(beta_k = perturb_beta_k(1.5),
                                   log_sigma_u = log(0.8)))
  )
  if (!isTRUE(multi_start)) start_sets <- start_sets[1]

  # NOTE: this best_kept/best_any + log-k blow-up guard parallels
  # .dd_tmb_multi_start() in R/dd-tmb.R. The choice family keeps its own copy
  # (choice-specific tmb_data, map, and beta_k/log_sigma_u perturbations, and to
  # avoid modifying the IP helper) — keep the selection logic in sync if either
  # changes.
  # R7: best_kept (passes log-k blow-up guard) AND best_any (lowest finite nll).
  best_kept <- NULL
  best_kept_nll <- Inf
  best_any <- NULL
  best_any_nll <- Inf
  opt_warnings <- character(0)
  for (s in start_sets) {
    res <- tryCatch({
      o <- TMB::MakeADFun(tmb_data, s, map = map, random = "u",
                          DLL = "beezdiscounting", silent = verbose < 2)
      opt_res <- .dd_tmb_run_optimizer(o, o$par, tmb_control, user_specified, verbose)
      list(obj = o, opt = opt_res$opt, nll = opt_res$opt$objective,
           warnings = opt_res$warnings)
    }, error = function(e) NULL)
    if (is.null(res) || !is.finite(res$nll)) next
    if (length(res$warnings)) opt_warnings <- c(opt_warnings, res$warnings)
    if (res$nll < best_any_nll) {
      best_any_nll <- res$nll
      best_any <- res
    }
    if (!.dd_logk_blowup(res$opt, tmb_data$X) && res$nll < best_kept_nll) {
      best_kept_nll <- res$nll
      best_kept <- res
    }
  }
  best <- best_kept
  if (is.null(best)) {
    if (is.null(best_any)) {
      stop("All starting value sets failed for fit_dd_choice().", call. = FALSE)
    }
    best <- best_any
    warning("All choice fits hit the log-k blow-up guard; returning the best ",
            "available fit (estimates may be unstable).", call. = FALSE)
  }
  obj <- best$obj
  opt <- best$opt
  converged <- isTRUE(opt$convergence == 0)

  try(obj$fn(opt$par), silent = TRUE)   # R4: refresh last.par.best so sdreport is fresh

  est <- .dd_choice_extract_estimates(obj, opt, prepared$n_subjects, intercept,
                                      verbose)
  # `family`/`equation` are reserved-but-cosmetic here: subject k is computed only
  # from beta_k + sigma_u * u_i (the choice likelihood is irrelevant to this helper).
  subject_pars <- .dd_tmb_compute_subject_pars(
    coefficients = est$coefficients, u_hat = est$u_hat,
    subject_levels = prepared$subject_levels, design_X = design$X,
    subject_id = prepared$subject_id, equation = equation, family = "sltb")
  # Contract: subject_pars is exactly id/u_i/k (subset/rename if the helper adds more).
  subject_pars <- subject_pars[, c("id", "u_i", "k")]

  nll <- opt$objective
  loglik <- -nll
  np <- length(opt$par)
  structure(list(
    call = cl, opt = opt,
    model = list(coefficients = est$coefficients, se = est$se),
    sdr = est$sdr, hessian_pd = est$hessian_pd,
    param_info = list(
      mode = "structural", equation = equation, intercept = isTRUE(intercept),
      n_obs = prepared$n_obs, n_subjects = prepared$n_subjects,
      n_random_effects = 1L, subject_levels = prepared$subject_levels,
      id_var = "id", x_var = "delay", y_var = "choice",
      factors = intersect(factors, all.vars(design$rhs)),
      factor_interaction = factor_interaction,
      continuous_covariates = continuous_covariates),
    formula_details = list(X = design$X, rhs = design$rhs,
                           contrasts = design$contrasts),
    subject_pars = subject_pars, loglik = loglik,
    AIC = 2 * nll + 2 * np, BIC = 2 * nll + np * log(prepared$n_obs),
    converged = converged, opt_warnings = opt_warnings,
    se_available = !is.null(est$sdr) && isTRUE(est$hessian_pd),
    data = prepared$data, coercion_info = validated$coercion_info
  ), class = "beezdiscounting_choice")
}
