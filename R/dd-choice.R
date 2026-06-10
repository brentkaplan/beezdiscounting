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
#' @param data Canonical frame (+ extras) from `.dd_validate_choice()`.
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
# Descriptive choice helpers (Young 2018) — pure R, no TMB.
# ==============================================================================

#' Default descriptive (Young 2018) predictor formula
#' @keywords internal
#' @noRd
.dd_choice_default_predictors <- function() {
  stats::as.formula("~ 0 + log(ll_amount / ss_amount) + log(delay + 1)")
}

#' Build the descriptive fixed (`Z`) and random-slope (`Zre`) designs
#'
#' `Z` comes from `predictors` (default = Young's two scale-invariant predictors,
#' no intercept) and is rank-guarded. `Zre` is ALWAYS Young's two predictors
#' (q = 2) when `random_slopes = TRUE`, else a 0-column matrix (q = 0, pooled).
#' Random slopes live on exactly the two Young predictors this slice
#' (general-q / re_formula deferred).
#'
#' @param data Cleaned canonical frame (`id`/`ss_amount`/`ll_amount`/`delay`/
#'   `choice` + any extras).
#' @param predictors One-sided formula for the fixed design, or `NULL` (Young).
#' @param random_slopes Logical; `TRUE` => q = 2 correlated slopes.
#' @return list(Z, Zre, q, predictors = formula, re_terms = colnames(Zre),
#'   contrasts, xlevels) — `contrasts`/`xlevels` let `predict()` rebuild a
#'   column-aligned `Z` on newdata and reject unseen factor levels.
#' @keywords internal
#' @noRd
.dd_choice_descriptive_design <- function(data, predictors = NULL,
                                          random_slopes = TRUE) {
  if (is.null(predictors)) predictors <- .dd_choice_default_predictors()
  if (!inherits(predictors, "formula")) {
    cli::cli_abort("{.arg predictors} must be a one-sided formula or {.code NULL}.")
  }
  mf <- stats::model.frame(predictors, data = data)
  Z <- stats::model.matrix(predictors, data = mf)
  if (ncol(Z) > 0L) {
    qrz <- qr(Z)
    if (qrz$rank < ncol(Z)) {
      aliased <- colnames(Z)[qrz$pivot[(qrz$rank + 1L):ncol(Z)]]
      cli::cli_abort(c(
        "The descriptive fixed-effect design is rank-deficient \\
         ({qrz$rank}/{ncol(Z)} columns independent).",
        "x" = "Non-estimable (aliased/collinear) column{?s}: {.val {aliased}}.",
        "i" = "Drop the redundant predictor or collapse empty cells."
      ))
    }
  }
  re_form <- .dd_choice_default_predictors()       # Young's 2, always
  Zre_full <- stats::model.matrix(re_form, data = data)
  if (isTRUE(random_slopes)) {
    Zre <- Zre_full
    q <- ncol(Zre)                                  # == 2L by construction
    # Random slopes need a matching fixed effect in Z (else the population-mean
    # slope is constrained to 0). The default `predictors` includes Young's two
    # terms; warn if a caller override drops them.
    missing_anchor <- setdiff(colnames(Zre), colnames(Z))
    if (length(missing_anchor)) {
      cli::cli_warn(c(
        "Random slope term(s) have no matching fixed effect in {.arg predictors}: \\
         {.val {missing_anchor}}.",
        "i" = "Their population-mean slope is constrained to 0. Include these \\
               terms in {.arg predictors} (the default does) unless intended."))
    }
  } else {
    Zre <- Zre_full[, integer(0), drop = FALSE]
    q <- 0L
  }
  list(Z = Z, Zre = Zre, q = q, predictors = predictors,
       re_terms = colnames(Zre),
       contrasts = attr(Z, "contrasts"),
       xlevels = stats::.getXlevels(stats::terms(predictors), mf))
}

#' 2-RE covariance from log-SDs + an unconstrained correlation
#'
#' Mirrors the `beezdemand` 2-RE pattern: `rho = tanh(cor_re)`,
#' `Sigma = diag(sd) %*% R %*% diag(sd)`, `L = chol(Sigma)` (lower-triangular).
#' For q = 2 `cor_re` is a scalar.
#'
#' @param log_sd_re Length-2 numeric (log SDs).
#' @param cor_re Length-1 numeric (unconstrained; `rho = tanh(cor_re)`).
#' @return list(Sigma = 2x2, L = lower-tri Cholesky, sd = c(sd1, sd2), rho).
#' @keywords internal
#' @noRd
.dd_chol_sigma <- function(log_sd_re, cor_re) {
  sd <- exp(log_sd_re)
  rho <- tanh(cor_re[[1]])
  Sigma <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                    rho * sd[1] * sd[2], sd[2]^2), nrow = 2L)
  L <- t(chol(Sigma))                                # lower-triangular
  list(Sigma = Sigma, L = L, sd = sd, rho = rho)
}

#' Descriptive (Young) linear predictor reference (pure R; compile-gate ref)
#'
#' `eta_i = Z_i theta + Zre_i (L b_{subj_i})`. `b` is the n_subjects x q matrix
#' of STANDARDIZED deviates; `L` is the Cholesky factor of Sigma.
#'
#' @keywords internal
#' @noRd
.dd_choice_descriptive_eta <- function(Z, theta, Zre, b, L, subject_id) {
  eta <- as.numeric(Z %*% theta)
  if (ncol(Zre) > 0L) {
    re <- t(L %*% t(b))                              # n_subjects x q on natural scale
    eta <- eta + rowSums(Zre * re[subject_id + 1L, , drop = FALSE])
  }
  unname(eta)
}


# ==============================================================================
# Structural choice fit pipeline (direct-k binomial GLMM via TMB).
# ==============================================================================

#' Full per-mode TMB `map` for the choice template
#'
#' Fixes every parameter block irrelevant to the active mode (thread through
#' EVERY `MakeADFun` call). Structural fixes the descriptive blocks
#' (`theta`/`log_sd_re`/`cor_re`/`b`) and `beta0` when intercept is off.
#' Descriptive fixes the structural blocks (`beta_k`/`log_sigma_u`/`log_gamma`/
#' `beta0`/`u`); with `random_slopes = FALSE` it additionally fixes the RE blocks
#' (`log_sd_re`/`cor_re`/`b`).
#' @keywords internal
#' @noRd
.dd_choice_build_map <- function(mode, intercept = FALSE,
                                 random_slopes = TRUE, starts) {
  na_f <- function(x) factor(rep(NA, length(unlist(x))))
  map <- list()
  if (mode == "structural") {
    if (!isTRUE(intercept)) map$beta0 <- factor(NA)
    map$theta <- na_f(starts$theta)
    map$log_sd_re <- na_f(starts$log_sd_re)
    map$cor_re <- na_f(starts$cor_re)
    map$b <- na_f(starts$b)
  } else {
    map$beta_k <- na_f(starts$beta_k)
    map$log_sigma_u <- factor(NA)
    map$log_gamma <- factor(NA)
    map$beta0 <- factor(NA)
    map$u <- na_f(starts$u)
    if (!isTRUE(random_slopes)) {
      map$log_sd_re <- na_f(starts$log_sd_re)
      map$cor_re <- na_f(starts$cor_re)
      map$b <- na_f(starts$b)
    }
  }
  map
}

#' Default starts for the structural choice model
#'
#' Includes the (unused-in-structural) descriptive blocks so the shared template
#' — which declares them unconditionally — can be instantiated; they are fixed
#' via [.dd_choice_build_map()].
#' @keywords internal
#' @noRd
.dd_choice_default_starts <- function(prepared, design) {
  p <- ncol(design$X)
  beta_k <- rep(0, p)
  beta_k[1] <- log(0.01)                  # generic mid-range k start
  list(
    beta_k = beta_k, log_sigma_u = log(0.5),
    log_gamma = log(1), beta0 = 0,
    u = matrix(0, nrow = prepared$n_subjects, ncol = 1L),
    theta = 0, log_sd_re = rep(log(0.5), 2L), cor_re = 0,
    b = matrix(0, nrow = prepared$n_subjects, ncol = 2L)
  )
}

#' Default start list carrying BOTH modes' parameter blocks (descriptive fit)
#' @keywords internal
#' @noRd
.dd_choice_full_starts <- function(prepared, n_x, n_z, q) {
  beta_k <- rep(0, max(n_x, 1L)); beta_k[1] <- log(0.01)
  list(
    beta_k = beta_k, log_sigma_u = log(0.5), log_gamma = log(1), beta0 = 0,
    u = matrix(0, prepared$n_subjects, 1L),
    theta = rep(0, max(n_z, 1L)),
    log_sd_re = rep(log(0.5), 2L),       # q = 2 slots; mapped out when q != 2
    cor_re = 0,
    b = matrix(0, prepared$n_subjects, 2L)
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
  if (any(c("theta", "log_sd_re", "cor_re") %in% par_names)) {
    stop("internal: descriptive parameters are free in a structural fit; check map threading.",
         call. = FALSE)
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

#' Extract estimates from a descriptive choice fit (sdreport + SEs + slopes)
#' @keywords internal
#' @noRd
.dd_choice_extract_descriptive <- function(obj, opt, prepared, design,
                                           random_slopes, verbose = 1) {
  sdr <- tryCatch(TMB::sdreport(obj), error = function(e) NULL)
  hessian_pd <- if (!is.null(sdr)) isTRUE(sdr$pdHess) else NA
  if (isTRUE(verbose >= 1) && (is.null(sdr) || !isTRUE(hessian_pd))) {
    cli::cli_warn(c("Standard errors may be unreliable (sdreport failed or \\
                     non-PD Hessian).", "i" = "Fixed-effect SEs/CIs will be {.val NA}."))
  }
  par_full <- opt$par
  par_names <- names(par_full)
  se_vec <- stats::setNames(rep(NA_real_, length(par_full)), par_names)
  if (!is.null(sdr)) {
    fs <- summary(sdr, "fixed")
    for (nm in unique(par_names)) {
      idx <- which(par_names == nm)
      rows <- fs[rownames(fs) == nm, , drop = FALSE]
      if (length(idx) && nrow(rows) == length(idx)) se_vec[idx] <- rows[, "Std. Error"]
    }
  }
  # Per-subject slopes b_i (natural scale) = L %*% b_std_i. Pooled (q=0): an
  # id-only frame so subject_pars is ALWAYS a data frame (never NULL).
  slopes <- data.frame(id = prepared$subject_levels, stringsAsFactors = FALSE)
  Sigma <- NULL
  if (isTRUE(random_slopes)) {
    cs <- .dd_chol_sigma(unname(par_full[par_names == "log_sd_re"]),
                         unname(par_full[par_names == "cor_re"]))
    Sigma <- cs$Sigma
    re_summary <- tryCatch(summary(sdr, "random"), error = function(e) NULL)
    b_std <- if (!is.null(re_summary)) {
      matrix(re_summary[, "Estimate"], nrow = prepared$n_subjects, ncol = 2L)
    } else {
      matrix(0, prepared$n_subjects, 2L)
    }
    nat <- t(cs$L %*% t(b_std))
    slopes <- data.frame(
      id = prepared$subject_levels,
      b_mag = nat[, 1], b_delay = nat[, 2], stringsAsFactors = FALSE)
  }
  list(coefficients = par_full, se = se_vec, sdr = sdr,
       hessian_pd = hessian_pd, Sigma = Sigma, slopes = slopes)
}

#' Descriptive sanity guard: reject non-finite params or runaway |theta|
#' @keywords internal
#' @noRd
.dd_choice_descr_blowup <- function(opt, max_abs = 30) {
  par <- opt$par
  v <- unname(par[names(par) %in% c("theta", "log_sd_re", "cor_re")])
  if (any(!is.finite(v))) return(TRUE)
  th <- unname(par[names(par) == "theta"])
  any(abs(th) > max_abs)
}

#' Descriptive (Young 2018) choice fit (internal dispatch target)
#' @keywords internal
#' @noRd
.fit_dd_choice_descriptive <- function(data, cl, id_var, ss_var, ll_var,
                                       delay_var, choice_var, predictors,
                                       random_slopes, factors,
                                       continuous_covariates, start_values,
                                       tmb_control, multi_start, verbose) {
  if (!is.null(factors) || !is.null(continuous_covariates)) {
    cli::cli_warn(c(
      "{.arg factors} / {.arg continuous_covariates} are ignored for \\
       {.code mode = \"descriptive\"}.",
      "i" = "The descriptive (Young 2018) model uses {.arg predictors} for its \\
             fixed design. Add between-subject terms there, e.g. \\
             {.code predictors = ~ 0 + log(ll_amount/ss_amount) + log(delay + 1) + group}."))
  }
  extra_cols <- unique(c(factors, continuous_covariates,
                         all.vars(predictors %||% .dd_choice_default_predictors())))
  extra_cols <- setdiff(extra_cols, c("ss_amount", "ll_amount", "delay"))
  validated <- .dd_validate_choice(data, id_var = id_var, ss_var = ss_var,
                                   ll_var = ll_var, delay_var = delay_var,
                                   choice_var = choice_var, extra_cols = extra_cols)
  prepared <- .dd_choice_prepare_data(validated$data, extra_cols = extra_cols)
  design <- .dd_choice_descriptive_design(prepared$data, predictors = predictors,
                                          random_slopes = random_slopes)

  tmb_data <- list(
    model = "ChoiceDiscounting", mode = 1L, eqn_type = 0L, has_intercept = 0L,
    choice = as.numeric(prepared$choice),
    ss_amount = as.numeric(prepared$ss_amount),
    ll_amount = as.numeric(prepared$ll_amount),
    delay = as.numeric(prepared$delay),
    subject_id = as.integer(prepared$subject_id),
    X = matrix(1, prepared$n_obs, 1L),     # structural design unused (mapped out)
    Z = as.matrix(design$Z), Zre = as.matrix(design$Zre),
    n_obs = as.integer(prepared$n_obs),
    n_subjects = as.integer(prepared$n_subjects),
    n_re = as.integer(design$q))

  starts <- .dd_choice_full_starts(prepared, n_x = 1L, n_z = ncol(design$Z),
                                   q = design$q)
  if (!is.null(start_values)) {
    for (nm in names(start_values)) if (nm %in% names(starts)) starts[[nm]] <- start_values[[nm]]
  }
  map <- .dd_choice_build_map("descriptive", random_slopes = random_slopes,
                              starts = starts)
  random_arg <- if (isTRUE(random_slopes)) "b" else NULL

  default_control <- list(iter_max = 1000, eval_max = 2000, optimizer = "nlminb",
                          rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
  user_specified <- names(tmb_control)
  tmb_control <- utils::modifyList(default_control, tmb_control)

  perturb <- function(d_theta, d_sd) utils::modifyList(
    starts, list(theta = starts$theta + d_theta, log_sd_re = log(d_sd)))
  start_sets <- list(starts, perturb(0.5, c(0.3, 0.2)), perturb(-0.5, c(0.8, 0.6)))
  if (!isTRUE(multi_start)) start_sets <- start_sets[1]

  best_kept <- NULL; best_kept_nll <- Inf
  best_any <- NULL; best_any_nll <- Inf
  opt_warnings <- character(0)
  for (s in start_sets) {
    res <- tryCatch({
      o <- TMB::MakeADFun(tmb_data, s, map = map, random = random_arg,
                          DLL = "beezdiscounting", silent = verbose < 2)
      opt_res <- .dd_tmb_run_optimizer(o, o$par, tmb_control, user_specified, verbose)
      list(obj = o, opt = opt_res$opt, nll = opt_res$opt$objective,
           warnings = opt_res$warnings)
    }, error = function(e) NULL)
    if (is.null(res) || !is.finite(res$nll)) next
    if (length(res$warnings)) opt_warnings <- c(opt_warnings, res$warnings)
    if (res$nll < best_any_nll) { best_any_nll <- res$nll; best_any <- res }
    if (!.dd_choice_descr_blowup(res$opt) && res$nll < best_kept_nll) {
      best_kept_nll <- res$nll; best_kept <- res
    }
  }
  best <- best_kept %||% best_any
  if (is.null(best)) stop("All starting value sets failed for the descriptive fit_dd_choice().",
                          call. = FALSE)
  obj <- best$obj; opt <- best$opt
  converged <- isTRUE(opt$convergence == 0)
  try(obj$fn(opt$par), silent = TRUE)

  free <- names(opt$par)
  if (any(c("beta_k", "log_gamma", "beta0", "log_sigma_u") %in% free)) {
    stop("internal: structural parameters are free in a descriptive fit; check map threading.",
         call. = FALSE)
  }
  est <- .dd_choice_extract_descriptive(obj, opt, prepared, design,
                                        random_slopes, verbose)
  nll <- opt$objective; np <- length(opt$par)
  structure(list(
    call = cl, opt = opt,
    model = list(coefficients = est$coefficients, se = est$se),
    sdr = est$sdr, hessian_pd = est$hessian_pd,
    param_info = list(
      mode = "descriptive", equation = NA_character_, intercept = FALSE,
      n_obs = prepared$n_obs, n_subjects = prepared$n_subjects,
      n_random_effects = design$q, subject_levels = prepared$subject_levels,
      id_var = "id", x_var = "delay", y_var = "choice",
      predictors = design$predictors, random_slopes = isTRUE(random_slopes),
      re_terms = design$re_terms, factors = NULL, factor_interaction = FALSE,
      continuous_covariates = continuous_covariates),
    formula_details = list(Z = design$Z, Zre = design$Zre,
                           predictors = design$predictors,
                           contrasts = design$contrasts, xlevels = design$xlevels),
    Sigma = est$Sigma, subject_pars = est$slopes,
    loglik = -nll, AIC = 2 * nll + 2 * np, BIC = 2 * nll + np * log(prepared$n_obs),
    converged = converged, opt_warnings = opt_warnings,
    se_available = !is.null(est$sdr) && isTRUE(est$hessian_pd),
    data = prepared$data, coercion_info = validated$coercion_info
  ), class = "beezdiscounting_choice")
}

#' Fit a trial-level SS-vs-LL choice model (binomial GLMM) via TMB
#'
#' Fits a trial-level binary choice model from two complementary perspectives.
#' The structural model (Young 2018) estimates the discount rate `k` directly
#' from choices via the scale-invariant relative-value comparison
#' `logit P(LL) = beta0 + gamma * ((ll/ss) * D(k, delay) - 1)`,
#' `k = exp(X beta_k + sigma_u u)`. The descriptive model forgoes a discount
#' function and instead characterises choices via separate magnitude and delay
#' sensitivities with optional correlated per-subject random slopes.
#'
#' @details
#' The **structural** model (`mode = "structural"`) parameterises choices
#' through a classical discount function (Mazur hyperbolic or exponential) and
#' shares the IP family's `k`/emmeans contract: `k` estimates and estimated
#' marginal means are accessible via [get_dd_param_emms()] and
#' [get_dd_comparisons()].  The **descriptive** model (`mode = "descriptive"`)
#' follows Young (2018): it regresses binary choice on
#' `log(ll_amount / ss_amount)` (magnitude sensitivity) and `log(delay + 1)`
#' (delay sensitivity) with no assumed discount function.  When
#' `random_slopes = TRUE` (default) each subject receives correlated random
#' slopes on these two predictors; the primary inferential targets are the
#' random-effect covariance (via [nlme::VarCorr()]) and per-subject slopes
#' (via [nlme::ranef()]), not emmeans.
#'
#' @param data Trial-level data frame (see the `*_var` args).
#' @param mode `"structural"` (default) estimates the discount rate `k`
#'   directly from choices via a discount function (Mazur or exponential);
#'   shares the IP family's `k`/emmeans contract.  `"descriptive"` fits the
#'   Young (2018) correlated random-slope logistic model with separate
#'   magnitude and delay sensitivity predictors; the inferential targets are
#'   `VarCorr` / `ranef`, not emmeans.
#' @param id_var,ss_var,ll_var,delay_var,choice_var Column names.
#' @param equation `"mazur"` or `"exponential"` (structural mode only).
#' @param intercept Logical; include the choice-bias `beta0` (default `FALSE`).
#'   Structural only (the descriptive design carries its own intercept policy via
#'   `predictors`).
#' @param predictors Descriptive (`mode = "descriptive"`) only: a one-sided
#'   formula for the fixed-effect design on the logit of choosing the larger-later
#'   reward, or `NULL` for Young's (2018) two scale-invariant predictors
#'   (`~ 0 + log(ll_amount / ss_amount) + log(delay + 1)`). Ignored when
#'   `mode = "structural"`.
#' @param random_slopes Descriptive only: logical; `TRUE` (default) fits two
#'   correlated per-subject random slopes on Young's two predictors
#'   (`log(ll/ss)` and `log(delay + 1)`), giving a full bivariate random-effect
#'   covariance; `FALSE` fits a pooled fixed-effect logistic model with no
#'   random effects. Ignored when `mode = "structural"`.
#' @param factors,factor_interaction,continuous_covariates Between-subject design
#'   on `log k` (same semantics as [fit_dd_tmb()]). Structural only.
#' @param start_values,tmb_control,multi_start,verbose,... As in [fit_dd_tmb()].
#' @return An object of class `beezdiscounting_choice`.
#' @seealso [simulate_dd_choice()] for data generation; [nlme::VarCorr()] and
#'   [nlme::ranef()] for descriptive-mode random-effect output;
#'   [get_dd_param_emms()] and [get_dd_comparisons()] for structural-mode
#'   emmeans.
#' @examples
#' \donttest{
#' # Structural model: estimate a discount rate k from binary choices
#' sim_s <- simulate_dd_choice(n_subjects = 30, mode = "structural", seed = 1)
#' fit_s <- fit_dd_choice(sim_s, mode = "structural", equation = "mazur")
#' summary(fit_s)
#'
#' # Descriptive (Young 2018) model: correlated per-subject magnitude/delay slopes
#' sim_d <- simulate_dd_choice(n_subjects = 30, mode = "descriptive", seed = 1)
#' fit_d <- fit_dd_choice(sim_d, mode = "descriptive")
#' summary(fit_d)
#' nlme::VarCorr(fit_d)
#' }
#' @export
fit_dd_choice <- function(data, mode = c("structural", "descriptive"),
                          id_var = "id", ss_var = "ss_amount",
                          ll_var = "ll_amount", delay_var = "delay",
                          choice_var = "choice",
                          equation = c("mazur", "exponential"),
                          intercept = FALSE,
                          predictors = NULL, random_slopes = TRUE,
                          factors = NULL, factor_interaction = FALSE,
                          continuous_covariates = NULL,
                          start_values = NULL,
                          tmb_control = list(iter_max = 1000, eval_max = 2000),
                          multi_start = TRUE, verbose = 1, ...) {
  cl <- match.call()
  mode <- match.arg(mode)
  equation <- match.arg(equation)
  if (mode == "descriptive") {
    return(.fit_dd_choice_descriptive(
      data = data, cl = cl, id_var = id_var, ss_var = ss_var, ll_var = ll_var,
      delay_var = delay_var, choice_var = choice_var,
      predictors = predictors, random_slopes = random_slopes,
      factors = factors, continuous_covariates = continuous_covariates,
      start_values = start_values, tmb_control = tmb_control,
      multi_start = multi_start, verbose = verbose))
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
    # Descriptive (mode 1) design blocks are declared unconditionally by the
    # shared template; supply inert placeholders for the structural path.
    Z = matrix(0, nrow = prepared$n_obs, ncol = 1L),
    Zre = matrix(0, nrow = prepared$n_obs, ncol = 1L),
    n_re = 0L,
    n_obs = as.integer(prepared$n_obs),
    n_subjects = as.integer(prepared$n_subjects)
  )
  starts <- .dd_choice_default_starts(prepared, design)
  if (!is.null(start_values)) {
    for (nm in names(start_values)) if (nm %in% names(starts)) starts[[nm]] <- start_values[[nm]]
  }
  map <- .dd_choice_build_map("structural", intercept = intercept,
                              random_slopes = FALSE, starts = starts)

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
