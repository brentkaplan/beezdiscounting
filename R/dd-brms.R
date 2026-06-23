# brms tier: indifference-point fitter ----------------------------------------------
#
# fit_dd_brms(): Bayesian mixed-effects delay discounting via brms/Stan,
# mirroring fit_dd_tmb()'s data handling and vocabulary. Design spec:
# beezdemand internal_docs/design/DESIGN-brms-tier.md (sections 1.3, 2.6).
# Helper machinery is the beezdemand brms tier ported with the .dd_ prefix.

#' Posterior draws matrix
#' @noRd
.dd_brms_draws_matrix <- function(object) {
  posterior::as_draws_matrix(object$brmsfit)
}

#' Canonicalize design column names the way brms sanitizes draw names
#'
#' brms strips parentheses/spaces (and other specials) from design column
#' names when naming draws: "(Intercept)" -> "b_logk_Intercept". Colons in
#' interaction columns are KEPT. Both sides of any X-vs-draws comparison go
#' through this.
#' @noRd
.dd_brms_canon_coefname <- function(x) {
  gsub("[^[:alnum:].:_]", "", x)
}

#' Ordered b_logk draw variables, validated against the fitted design
#'
#' THE single alignment point between posterior draw names and the stored
#' TMB-style design matrix (Codex 048-B1: a raw grep() is positional and a
#' length check cannot catch reordering or name mangling). Resolution order:
#' the fit-time map `formula_details$logk_draw_vars` (captured from the brms
#' standata design, see `.dd_brms_logk_standata_map()`), else the canonical
#' name reconstruction from `colnames(formula_details$X)` for objects fit
#' before the map existed. Either way the result is validated for set
#' equality against the draw columns actually present and returned in
#' FITTED-DESIGN ORDER, so `draws[, vars]` is column-aligned with
#' `formula_details$X` by NAME, never by position.
#' @noRd
.dd_brms_logk_draw_vars <- function(object, draws_cols) {
  fitted_cols <- colnames(object$formula_details$X)
  expected <- object$formula_details$logk_draw_vars
  if (is.null(expected)) {
    expected <- paste0("b_logk_", .dd_brms_canon_coefname(fitted_cols))
  }
  found <- grep("^b_logk_", draws_cols, value = TRUE)
  if (
    length(expected) != length(fitted_cols) ||
      anyDuplicated(expected) > 0 ||
      !setequal(expected, found)
  ) {
    stop(
      "Internal error: cannot align b_logk_* posterior draws with the ",
      "fitted log-k design columns.\n  Design: ",
      paste(fitted_cols, collapse = ", "),
      "\n  Draws:  ",
      paste(found, collapse = ", "),
      call. = FALSE
    )
  }
  expected
}

#' Capture the authoritative b_logk draw-name map from the brms standata
#'
#' The Stan coefficient vector `b_logk` is ordered by the columns of the
#' standata design `X_logk`; brms names the draws
#' `paste0("b_logk_", colnames(X_logk))`. Cross-checks that order
#' elementwise (canonicalized) against the stored TMB-style design so any
#' brms-side reordering aborts AT FIT TIME instead of silently relabeling
#' coefficients downstream.
#' @noRd
.dd_brms_logk_standata_map <- function(brmsfit, X) {
  sdata_cols <- colnames(brms::standata(brmsfit)$X_logk)
  if (
    !identical(
      .dd_brms_canon_coefname(colnames(X)),
      .dd_brms_canon_coefname(sdata_cols)
    )
  ) {
    stop(
      "Internal error: the brms design for log k does not match the fitted ",
      "design matrix.\n  Fitted: ",
      paste(colnames(X), collapse = ", "),
      "\n  brms:   ",
      paste(sdata_cols, collapse = ", "),
      call. = FALSE
    )
  }
  paste0("b_logk_", sdata_cols)
}

#' Evaluate code under a temporary RNG seed, restoring the caller's stream
#' @noRd
.dd_brms_with_seed <- function(seed, code) {
  if (is.na(seed)) {
    return(code)
  }
  old_seed <- if (
    exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  ) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed + 1L)
  code # promise evaluates here, after set.seed()
}

#' MCMC diagnostics (all-NA-safe: chains = 1 must not fake convergence)
#' @noRd
.dd_brms_mcmc_diagnostics <- function(brmsfit, max_treedepth = 10) {
  sm <- posterior::summarise_draws(
    posterior::as_draws_array(brmsfit),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  np <- brms::nuts_params(brmsfit)
  n_div <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
  td <- np$Value[np$Parameter == "treedepth__"]
  safe_extreme <- function(x, fn) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else fn(x)
  }
  list(
    rhat_max = safe_extreme(sm$rhat, max),
    ess_bulk_min = safe_extreme(sm$ess_bulk, min),
    ess_tail_min = safe_extreme(sm$ess_tail, min),
    num_divergent = as.integer(n_div),
    num_max_treedepth = as.integer(sum(td >= max_treedepth, na.rm = TRUE))
  )
}

#' Convergence decision: all criteria finite AND passing
#' @noRd
.dd_brms_converged <- function(diag) {
  isTRUE(
    is.finite(diag$rhat_max) &&
      diag$rhat_max < 1.01 &&
      diag$num_divergent == 0 &&
      is.finite(diag$ess_bulk_min) &&
      diag$ess_bulk_min >= 400
  )
}

#' Transform draws between report spaces (estimation scale = natural log)
#' @noRd
.dd_brms_transform_draws <- function(draws, to = "natural") {
  switch(
    to,
    natural = exp(draws),
    log10 = draws / log(10),
    log = draws,
    internal = draws,
    stop("Unsupported report space: ", to, call. = FALSE)
  )
}

#' Fit a Bayesian Mixed-Effects Discounting Model via brms
#'
#' Fits the four discounting equations of [fit_dd_tmb()] (`"mazur"`,
#' `"exponential"`, `"green-myerson"`, `"rachlin"`) as Bayesian nonlinear
#' mixed-effects models, with the discount rate estimated on the natural-log
#' scale (`logk`, subject random intercept; `k ~ 1`, the v1 scope) and the
#' shape exponent (`logs`) population-level for the two-parameter equations.
#'
#' `family = "beta"` (default) uses `Beta(link = "identity")` with the mean
#' squished into `(1e-6, 1 - 1e-6)` -- the closest brms analog of the TMB
#' SLT-beta (`family = "sltb"` has no brms equivalent and errors with this
#' pointer). Boundary observations (`y` exactly 0 or 1) are handled per
#' `boundary`: `"squeeze"` (default) applies the Smithson-Verkuilen
#' transform `y* = (y (N - 1) + 0.5) / N` to all responses (message reports
#' the boundary count); `"zoib"` switches to `zero_one_inflated_beta`
#' (statistically more honest but changes the estimand -- k then describes
#' interior responses only); `"error"` refuses to fit.
#' `family = "gaussian"` matches `fit_dd_tmb(family = "gaussian")`
#' wherever the TMB template's mu clamp into `[1e-6, 1 - 1e-6]` does not
#' bind (everywhere except extreme decay underflow).
#'
#' @param data Long-format data frame (one row per subject-delay).
#' @param y_var,x_var,id_var Column names (canonical defaults), as in
#'   [fit_dd_tmb()].
#' @param equation Discounting equation.
#' @param family `"beta"` or `"gaussian"`; `"sltb"` errors with guidance.
#' @param boundary Boundary handling for the beta family (see Details).
#' @param random_effects v1 supports `k ~ 1` only; `k + phi ~ 1` (supported
#'   by the TMB tier) errors with the planned brms route.
#' @param factors,factor_interaction,continuous_covariates Fixed-effect
#'   design on `logk`, as in [fit_dd_tmb()].
#' @param ll,response_scale Response coercion, as in [fit_dd_tmb()].
#' @param prior Optional `brmsprior`; user rows override the defaults.
#' @param autoscale_priors Anchor the `logk` prior to the median delay (see
#'   [default_dd_priors()]).
#' @param chains,iter,warmup,thin,cores,seed,backend,control,sample_prior
#'   MCMC settings passed to [brms::brm()].
#' @param init `"prior_center"` (default), `"tmb"` (a quiet
#'   [fit_dd_tmb()] pre-fit supplies the centers, with prior_center
#'   fallback on failure; the beta family maps to the TMB sltb pre-fit),
#'   `"random"`, or `"0"`; or a list/function passed through to
#'   [brms::brm()].
#' @param loo Compute and store [brms::loo()] at fit time.
#' @param file,file_refit Passed to [brms::brm()] for fit caching.
#' @param verbose 0 (silent), 1 (messages), 2 (full Stan output).
#' @param ... Passed through to [brms::brm()].
#'
#' @return An object of class `beezdiscounting_brms`: `model` (posterior
#'   medians/SDs on the estimation scale under TMB names `beta_k`/`log_s`,
#'   plus `variance_components`), `brmsfit`, `subject_pars`
#'   (`id`, `k`, `k_lower`, `k_upper`), `converged` (Rhat < 1.01, no
#'   divergences, bulk ESS >= 400), `mcmc_info`, `loo`, `data`,
#'   `param_info`, `priors`, `autoscale_info`.
#'
#' @seealso [fit_dd_tmb()]; [default_dd_priors()].
#' @export
fit_dd_brms <- function(
  data,
  y_var = "y",
  x_var = "x",
  id_var = "id",
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  family = c("beta", "gaussian"),
  boundary = c("squeeze", "zoib", "error"),
  random_effects = k ~ 1,
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  ll = NULL,
  response_scale = c("proportion", "percent", "amount"),
  prior = NULL,
  autoscale_priors = TRUE,
  chains = 4,
  iter = 2000,
  warmup = floor(iter / 2),
  thin = 1,
  cores = getOption("mc.cores", 1L),
  seed = NA,
  backend = getOption("brms.backend", "rstan"),
  control = list(adapt_delta = 0.95),
  init = c("prior_center", "tmb", "random", "0"),
  sample_prior = "no",
  loo = TRUE,
  file = NULL,
  file_refit = getOption("brms.file_refit", "on_change"),
  verbose = 1,
  ...
) {
  .dd_brms_check_installed()
  cl <- match.call()
  equation <- match.arg(equation)
  if (!identical(family, "sltb")) {
    family <- match.arg(family)
  }
  boundary <- match.arg(boundary)
  response_scale <- match.arg(response_scale)

  # v1 random-effects scope: k ~ 1 only (the TMB tier's k + phi ~ 1 is a
  # planned brms extension via lf(phi ~ 1 + (1|id))).
  re <- .dd_normalize_re(random_effects, data = data)
  re_params <- re$blocks[[1]]$param
  if (!identical(re_params, "k")) {
    stop(
      "v1 of the brms tier supports `random_effects = k ~ 1` only.\n",
      "phi random effects (`k + phi ~ 1`) are planned via a distributional ",
      "formula once the sltb-analog question is settled; use fit_dd_tmb() ",
      "meanwhile.",
      call. = FALSE
    )
  }

  validated <- .dd_validate_ip(
    data,
    y_var = y_var,
    x_var = x_var,
    id_var = id_var,
    ll = ll,
    extra_cols = unique(c(factors, continuous_covariates)),
    response_scale = response_scale
  )
  prep <- .dd_tmb_prepare_data(
    validated$data,
    y_var = y_var,
    x_var = x_var,
    id_var = id_var,
    extra_cols = unique(c(factors, continuous_covariates))
  )
  d <- prep$data
  for (f in factors) {
    d[[f]] <- as.factor(d[[f]])
  }

  # The TMB design path carries the guards (rank-deficiency rejection,
  # covariate checks) the bare model.matrix() call would skip (Codex 039-R1).
  design <- .dd_tmb_build_design(
    d,
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates
  )

  spec <- .dd_brms_formula(
    equation = equation,
    family = family,
    boundary = boundary,
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    data = d
  )

  # Boundary handling (beta family only; gaussian models raw y).
  n_boundary <- sum(d$y <= 0 | d$y >= 1)
  if (family == "beta") {
    if (boundary == "error" && n_boundary > 0) {
      stop(
        n_boundary,
        " boundary response(s) (y = 0 or 1) present; ",
        "refusing to fit with boundary = \"error\". Use \"squeeze\" or \"zoib\".",
        call. = FALSE
      )
    }
    if (boundary == "squeeze") {
      n_y <- nrow(d)
      d$y <- (d$y * (n_y - 1) + 0.5) / n_y
      if (n_boundary > 0 && verbose >= 1) {
        message(
          "Smithson-Verkuilen squeeze applied to all ",
          n_y,
          " responses (",
          n_boundary,
          " were exactly 0 or 1)."
        )
      }
    }
  }

  if (equation == "rachlin") {
    d$xzero <- as.numeric(d$x == 0)
    d$xsafe <- ifelse(d$x == 0, 1, d$x)
  }

  defaults <- default_dd_priors(
    equation,
    family = family,
    data = if (isTRUE(autoscale_priors)) d else NULL,
    y_var = "y",
    x_var = "x",
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    autoscale = isTRUE(autoscale_priors)
  )
  merged_priors <- .dd_brms_merge_priors(prior, defaults)
  autoscale_info <- attr(defaults, "autoscale_info")

  brms::validate_prior(
    merged_priors,
    formula = spec$formula,
    data = d,
    family = spec$family
  )

  if (verbose >= 1 && !is.null(autoscale_info)) {
    message(
      "Default priors autoscaled to the data (median delay = ",
      signif(autoscale_info$median_delay, 4),
      "); see $autoscale_info and default_dd_priors(). ",
      "Disable with autoscale_priors = FALSE."
    )
  }

  init_obj <- .dd_brms_build_inits(
    init = init,
    spec = spec,
    data = d,
    chains = chains,
    seed = seed,
    autoscale_info = autoscale_info,
    family = family,
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    equation = equation
  )

  if (verbose >= 1) {
    message(
      "Compiling and sampling via ",
      backend,
      " (first call compiles the Stan model, typically 1-2 minutes; ",
      "pass file = to cache)."
    )
  }

  fit_call <- function() {
    brms::brm(
      formula = spec$formula,
      data = d,
      family = spec$family,
      prior = merged_priors,
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      cores = cores,
      seed = seed,
      backend = backend,
      control = control,
      init = init_obj,
      sample_prior = sample_prior,
      file = file,
      file_refit = file_refit,
      refresh = if (verbose >= 2) max(floor(iter / 10), 1) else 0,
      silent = if (verbose >= 2) 0 else 2,
      ...
    )
  }
  brmsfit <- if (verbose >= 1) fit_call() else suppressMessages(fit_call())

  .dd_brms_assemble_fit(
    brmsfit = brmsfit,
    spec = spec,
    data = d,
    design = design,
    coercion_info = validated$coercion_info,
    merged_priors = merged_priors,
    autoscale_info = autoscale_info,
    call = cl,
    mcmc_settings = list(
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      seed = seed,
      backend = backend,
      max_treedepth = if (!is.null(control$max_treedepth)) {
        control$max_treedepth
      } else {
        10
      }
    ),
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    family = family,
    boundary = boundary,
    response_scale = response_scale,
    n_boundary = n_boundary,
    compute_loo = isTRUE(loo),
    verbose = verbose
  )
}

#' Per-chain inits at prior centers (or a quiet TMB pre-fit)
#' @noRd
.dd_brms_build_inits <- function(
  init,
  spec,
  data,
  chains,
  seed,
  autoscale_info,
  family,
  factors,
  factor_interaction,
  continuous_covariates,
  equation
) {
  if (is.list(init) || is.function(init)) {
    return(init)
  }
  init <- match.arg(init, c("prior_center", "tmb", "random", "0"))
  if (identical(init, "random")) {
    return("random")
  }
  if (identical(init, "0")) {
    return(0)
  }

  logk_center <- if (!is.null(autoscale_info)) {
    -log(autoscale_info$median_delay)
  } else {
    -4.5
  }
  centers <- list(logk = logk_center, logs = 0, sd = 0.5, phi = 8, sigma = 0.1)

  if (identical(init, "tmb")) {
    tmb_centers <- tryCatch(
      {
        tmb_fit <- fit_dd_tmb(
          data,
          y_var = "y",
          x_var = "x",
          id_var = "id",
          equation = equation,
          family = if (identical(family, "beta")) "sltb" else "gaussian",
          factors = factors,
          factor_interaction = factor_interaction,
          continuous_covariates = continuous_covariates,
          verbose = 0
        )
        coefs <- tmb_fit$model$coefficients
        out <- list(
          # full fixed-effect vector: factor/covariate designs get every
          # coefficient centered at the TMB MLE, not just the intercept
          # (Codex 041-R2)
          beta_k_vec = unname(coefs[names(coefs) == "beta_k"]),
          logk = unname(coefs[names(coefs) == "beta_k"])[1],
          sd = exp(unname(coefs[["log_sigma_u"]]))
        )
        if ("log_s" %in% names(coefs)) {
          out$logs <- unname(coefs[["log_s"]])
        }
        if ("log_phi" %in% names(coefs)) {
          out$phi <- exp(unname(coefs[["log_phi"]]))
        }
        if ("log_sigma_e" %in% names(coefs)) {
          out$sigma <- exp(unname(coefs[["log_sigma_e"]]))
        }
        out
      },
      error = function(e) {
        warning(
          "TMB pre-fit for init = \"tmb\" failed (",
          conditionMessage(e),
          "); falling back to init = \"prior_center\".",
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(tmb_centers)) {
      centers[names(tmb_centers)] <- tmb_centers
    }
  }
  logk_center <- centers$logk

  rhs <- build_fixed_rhs(
    factors = factors,
    factor_interaction = factor_interaction,
    continuous_covariates = continuous_covariates,
    data = data
  )
  K <- ncol(stats::model.matrix(stats::as.formula(rhs), data = data))
  n_id <- length(unique(data$id))

  beta_k_center <- if (
    !is.null(centers$beta_k_vec) &&
      length(centers$beta_k_vec) == K
  ) {
    centers$beta_k_vec
  } else {
    c(logk_center, rep(0, K - 1))
  }

  jit <- function(n = 1) stats::rnorm(n, 0, 0.1)
  one_chain <- function() {
    out <- list(
      b_logk = as.array(beta_k_center + jit(K)),
      sd_1 = as.array(abs(centers$sd + jit())),
      z_1 = matrix(jit(n_id), nrow = 1)
    )
    if (spec$has_s) {
      out$b_logs <- as.array(centers$logs + jit())
    }
    if (identical(family, "beta")) {
      out$phi <- abs(centers$phi + stats::rnorm(1, 0, 1)) # TMB start: log_aux = log(8)
    } else {
      out$sigma <- abs(centers$sigma + jit())
    }
    out
  }
  .dd_brms_with_seed(seed, lapply(seq_len(chains), function(i) one_chain()))
}

#' Assemble the beezdiscounting_brms object
#' @noRd
.dd_brms_assemble_fit <- function(
  brmsfit,
  spec,
  data,
  design,
  coercion_info,
  merged_priors,
  autoscale_info,
  call,
  mcmc_settings,
  factors,
  factor_interaction,
  continuous_covariates,
  family,
  boundary,
  response_scale,
  n_boundary,
  compute_loo,
  verbose
) {
  X <- design$X

  obj <- structure(
    list(
      model = NULL,
      brmsfit = brmsfit,
      subject_pars = NULL,
      converged = NA,
      mcmc_info = NULL,
      loo = NULL,
      loglik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      data = data,
      coercion_info = coercion_info,
      param_info = list(
        equation = spec$equation,
        family = family,
        boundary = if (family == "beta") boundary else NA_character_,
        n_boundary = n_boundary,
        has_s = spec$has_s,
        has_phi = family == "beta",
        n_obs = nrow(data),
        n_subjects = length(unique(data$id)),
        n_random_effects = 1L,
        subject_levels = sort(unique(as.character(data$id))),
        id_var = "id",
        x_var = "x",
        y_var = "y",
        response_scale = response_scale,
        factors = factors,
        factor_interaction = factor_interaction,
        continuous_covariates = continuous_covariates
      ),
      formula_details = list(
        X = X,
        rhs = design$rhs,
        contrasts = design$contrasts,
        brmsformula = spec$formula,
        family = spec$family
      ),
      priors = merged_priors,
      autoscale_info = autoscale_info,
      call = call
    ),
    class = c("beezdiscounting_brms", "list")
  )

  # Authoritative draw-name map captured at fit time + ordered, validated
  # alignment of every b_logk_* draw with the design columns (Codex 048-B1).
  obj$formula_details$logk_draw_vars <- .dd_brms_logk_standata_map(brmsfit, X)
  draws <- .dd_brms_draws_matrix(obj)
  k_vars <- .dd_brms_logk_draw_vars(obj, colnames(draws))
  vars <- k_vars
  names_out <- rep("beta_k", length(k_vars))
  if (spec$has_s) {
    vars <- c(vars, "b_logs_Intercept")
    names_out <- c(names_out, "log_s")
  }
  est <- apply(draws[, vars, drop = FALSE], 2, stats::median)
  se <- apply(draws[, vars, drop = FALSE], 2, stats::sd)
  names(est) <- names_out
  names(se) <- names_out

  ln10 <- log(10)
  vc_rows <- list(data.frame(
    Component = "sigma_u (log10-k RE SD)",
    Estimate = stats::median(as.numeric(draws[, "sd_id__logk_Intercept"])) /
      ln10,
    Scale = "log10",
    stringsAsFactors = FALSE
  ))
  if (family == "beta") {
    vc_rows[[2]] <- data.frame(
      Component = "phi (precision)",
      Estimate = stats::median(as.numeric(draws[, "phi"])),
      Scale = "natural",
      stringsAsFactors = FALSE
    )
  } else {
    vc_rows[[2]] <- data.frame(
      Component = "sigma_e (Residual SD)",
      Estimate = stats::median(as.numeric(draws[, "sigma"])),
      Scale = "natural",
      stringsAsFactors = FALSE
    )
  }

  obj$model <- list(
    coefficients = est,
    se = se,
    variance_components = do.call(rbind, vc_rows)
  )

  obj$subject_pars <- .dd_brms_subject_pars(obj, draws = draws)

  diag <- .dd_brms_mcmc_diagnostics(
    brmsfit,
    max_treedepth = mcmc_settings$max_treedepth
  )
  obj$mcmc_info <- c(
    mcmc_settings[c("chains", "iter", "warmup", "thin", "seed", "backend")],
    list(brms_version = as.character(utils::packageVersion("brms"))),
    diag
  )
  obj$converged <- .dd_brms_converged(diag)

  if (isTRUE(compute_loo)) {
    obj$loo <- tryCatch(
      brms::loo(brmsfit),
      error = function(e) {
        warning("loo computation failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  }

  if (verbose >= 1 && !obj$converged) {
    warning(
      sprintf(
        "MCMC diagnostics flag potential non-convergence (rhat_max = %.4g, divergences = %d, min bulk ESS = %d).",
        diag$rhat_max,
        diag$num_divergent,
        round(diag$ess_bulk_min)
      ),
      call. = FALSE
    )
  }

  obj
}

#' Per-subject k summaries from posterior draws
#'
#' `k_i = exp(x_i' beta_k + r_i)` per draw, using each subject's first
#' design row (the v1 fixed-effect design is between-subject by
#' construction; within-subject designs return NULL with a warning).
#' @noRd
.dd_brms_subject_pars <- function(
  object,
  draws = NULL,
  probs = c(0.025, 0.975)
) {
  if (is.null(draws)) {
    draws <- .dd_brms_draws_matrix(object)
  }
  dat <- object$data
  ids <- object$param_info$subject_levels
  X_full <- object$formula_details$X

  first_row <- match(ids, as.character(dat$id))
  for (i in seq_along(ids)) {
    rows <- X_full[as.character(dat$id) == ids[i], , drop = FALSE]
    if (nrow(unique(rows)) > 1L) {
      warning(
        "Subject-level k undefined: the fixed-effect design varies within subject ",
        ids[i],
        "; subject_pars skipped.",
        call. = FALSE
      )
      return(NULL)
    }
  }

  k_vars <- .dd_brms_logk_draw_vars(object, colnames(draws))
  b <- as.matrix(draws[, k_vars, drop = FALSE])
  r_vars <- paste0("r_id__logk[", ids, ",Intercept]")
  missing <- setdiff(r_vars, colnames(draws))
  if (length(missing) > 0) {
    stop(
      "Internal error: random-effect draws not found: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  r <- as.matrix(draws[, r_vars, drop = FALSE])

  X_sub <- X_full[first_row, , drop = FALSE]
  logk_i <- b %*% t(X_sub) + r
  k_i <- exp(matrix(as.numeric(logk_i), nrow = nrow(logk_i)))

  data.frame(
    id = ids,
    k = apply(k_i, 2, stats::median),
    k_lower = apply(k_i, 2, stats::quantile, probs = probs[1]),
    k_upper = apply(k_i, 2, stats::quantile, probs = probs[2]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
