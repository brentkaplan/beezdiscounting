# brms tier: structural choice fitter ----------------------------------------------
#
# fit_dd_choice_brms(): Bayesian trial-level choice model mirroring
# fit_dd_choice(mode = "structural"). Design: beezdemand
# internal_docs/design/DESIGN-brms-tier.md sections 1.4, 2.7.

#' Fit a Bayesian Structural Choice Discounting Model via brms
#'
#' The TMB structural likelihood of [fit_dd_choice()]:
#' `logit P(LL) = [b0] + gamma * ((ll/ss) * D(k, delay) - 1)` with choice
#' sensitivity `gamma = exp(loggamma)` and the discount rate `k =
#' exp(logk)` carrying the subject random intercept. With
#' `bernoulli("logit")` the brms nonlinear formula IS the logit, so the
#' likelihood matches TMB exactly.
#'
#' v1 implements `mode = "structural"` only: the descriptive (Young 2018)
#' model is a plain logistic GLMM expressible directly with
#' `brms::brm()`.
#'
#' @param data Long-format trial-level data (one row per choice).
#' @param mode `"structural"` (only mode in v1).
#' @param id_var,ss_var,ll_var,delay_var,choice_var Column names, as in
#'   [fit_dd_choice()].
#' @param equation `"mazur"` or `"exponential"`.
#' @param intercept Include the logit-scale bias term `b0`.
#' @inheritParams fit_dd_brms
#' @return An object of class `beezdiscounting_choice_brms`. Coefficients
#'   are posterior medians on the estimation scale under the TMB names
#'   (`beta_k`, `log_gamma`, and `beta0` when `intercept = TRUE`).
#' @seealso [fit_dd_choice()]; [default_dd_choice_priors()].
#' @export
fit_dd_choice_brms <- function(
  data,
  mode = c("structural", "descriptive"),
  id_var = "id",
  ss_var = "ss_amount",
  ll_var = "ll_amount",
  delay_var = "delay",
  choice_var = "choice",
  equation = c("mazur", "exponential"),
  intercept = FALSE,
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
  init = c("prior_center", "random", "0"),
  sample_prior = "no",
  loo = TRUE,
  file = NULL,
  file_refit = getOption("brms.file_refit", "on_change"),
  verbose = 1,
  ...
) {
  .dd_brms_check_installed()
  cl <- match.call()
  mode <- match.arg(mode)
  equation <- match.arg(equation)

  if (mode == "descriptive") {
    stop(
      "mode = \"descriptive\" is not implemented in the brms tier: the ",
      "Young (2018) descriptive model is a plain logistic GLMM you can fit ",
      "directly with brms::brm() (or use fit_dd_choice()).",
      call. = FALSE
    )
  }

  validated <- .dd_validate_choice(
    data,
    id_var = id_var, ss_var = ss_var, ll_var = ll_var,
    delay_var = delay_var, choice_var = choice_var
  )
  prepared <- .dd_choice_prepare_data(validated$data)
  d <- prepared$data
  d$rel <- d$ll_amount / d$ss_amount
  d$id <- droplevels(as.factor(d$id))

  spec <- .dd_brms_choice_formula(equation = equation, intercept = intercept)

  defaults <- default_dd_choice_priors(
    equation,
    intercept = intercept,
    data = if (isTRUE(autoscale_priors)) d else NULL,
    delay_var = "delay",
    autoscale = isTRUE(autoscale_priors)
  )
  merged_priors <- .dd_brms_merge_priors(prior, defaults)
  autoscale_info <- attr(defaults, "autoscale_info")

  brms::validate_prior(
    merged_priors,
    formula = spec$formula, data = d, family = spec$family
  )

  init_obj <- .dd_brms_build_choice_inits(
    init = init,
    chains = chains,
    seed = seed,
    autoscale_info = autoscale_info,
    intercept = intercept,
    n_id = length(unique(d$id))
  )

  if (verbose >= 1) {
    message(
      "Compiling and sampling via ", backend,
      " (first call compiles the Stan model; pass file = to cache)."
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
      data = d,
      param_info = list(
        mode = "structural",
        equation = equation,
        intercept = intercept,
        n_obs = nrow(d),
        n_subjects = length(unique(d$id)),
        n_random_effects = 1L,
        subject_levels = sort(unique(as.character(d$id)))
      ),
      formula_details = list(
        X = matrix(1, nrow(d), 1, dimnames = list(NULL, "(Intercept)")),
        brmsformula = spec$formula,
        family = spec$family
      ),
      priors = merged_priors,
      autoscale_info = autoscale_info,
      call = cl
    ),
    class = c("beezdiscounting_choice_brms", "list")
  )

  draws <- .dd_brms_draws_matrix(obj)
  vars <- c("b_logk_Intercept", "b_loggamma_Intercept")
  names_out <- c("beta_k", "log_gamma")
  if (isTRUE(intercept)) {
    vars <- c(vars, "b_b0_Intercept")
    names_out <- c(names_out, "beta0")
  }
  est <- apply(draws[, vars, drop = FALSE], 2, stats::median)
  se <- apply(draws[, vars, drop = FALSE], 2, stats::sd)
  names(est) <- names_out
  names(se) <- names_out

  obj$model <- list(
    coefficients = est,
    se = se,
    variance_components = data.frame(
      Component = "sigma_u (log10-k RE SD)",
      Estimate = stats::median(as.numeric(draws[, "sd_id__logk_Intercept"])) / log(10),
      Scale = "log10",
      stringsAsFactors = FALSE
    )
  )
  obj$subject_pars <- .dd_brms_subject_pars(obj, draws = draws)

  diag <- .dd_brms_mcmc_diagnostics(
    brmsfit,
    max_treedepth = if (!is.null(control$max_treedepth)) control$max_treedepth else 10
  )
  obj$mcmc_info <- c(
    list(
      chains = chains, iter = iter, warmup = warmup, thin = thin,
      seed = seed, backend = backend,
      brms_version = as.character(utils::packageVersion("brms"))
    ),
    diag
  )
  obj$converged <- .dd_brms_converged(diag)

  if (isTRUE(loo)) {
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
        diag$rhat_max, diag$num_divergent, round(diag$ess_bulk_min)
      ),
      call. = FALSE
    )
  }

  obj
}

#' Per-chain inits for the choice model
#' @noRd
.dd_brms_build_choice_inits <- function(
  init, chains, seed, autoscale_info, intercept, n_id
) {
  if (is.list(init) || is.function(init)) {
    return(init)
  }
  init <- match.arg(init, c("prior_center", "random", "0"))
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
  jit <- function(n = 1) stats::rnorm(n, 0, 0.1)
  one_chain <- function() {
    out <- list(
      b_logk = as.array(logk_center + jit()),
      b_loggamma = as.array(1 + jit()), # prior center: loggamma ~ normal(1, 1)
      sd_1 = as.array(abs(0.5 + jit())),
      z_1 = matrix(jit(n_id), nrow = 1)
    )
    if (isTRUE(intercept)) {
      out$b_b0 <- as.array(jit())
    }
    out
  }
  .dd_brms_with_seed(seed, lapply(seq_len(chains), function(i) one_chain()))
}

# --- methods -------------------------------------------------------------------------

#' Coefficient table for the choice model (draws-based)
#'
#' Mirrors `tidy.beezdiscounting_choice()`: the log-k coefficient is the
#' `"fixed"` row; `gamma` (choice sensitivity) and `beta0` (choice bias)
#' are `"shape"` rows; `beta0` is on the identity (logit-intercept) scale
#' and is NEVER transformed across report spaces.
#' @noRd
.dd_brms_choice_coef_table <- function(object, report_space = "natural") {
  report_space <- match.arg(report_space, c("natural", "log10", "internal", "log"))
  to <- if (report_space == "internal") "log" else report_space
  draws <- .dd_brms_draws_matrix(object)

  row_for <- function(var, term, transform, scale_label, component,
                      display = NULL) {
    dj <- as.numeric(draws[, var])
    tdraws <- if (transform) .dd_brms_transform_draws(dj, to = to) else dj
    tibble::tibble(
      term = term,
      estimate = stats::median(tdraws),
      std.error = stats::sd(tdraws),
      statistic = NA_real_,
      p.value = NA_real_,
      component = component,
      estimate_scale = if (transform) to else scale_label,
      term_display = if (!is.null(display)) {
        display
      } else if (transform) {
        .dd_term_display_space(term, to)
      } else {
        term
      }
    )
  }

  g_disp <- switch(to,
    natural = "gamma",
    log10 = "log10(gamma)",
    log = "log(gamma)"
  )
  rows <- list(
    row_for("b_logk_Intercept", "k:(Intercept)", TRUE, NULL, "fixed"),
    row_for("b_loggamma_Intercept", "gamma", TRUE, NULL, "shape",
      display = g_disp
    )
  )
  if (isTRUE(object$param_info$intercept)) {
    rows <- c(rows, list(
      row_for("b_b0_Intercept", "beta0", FALSE, "identity", "shape")
    ))
  }
  dplyr::bind_rows(rows)
}

#' Credible intervals for a beezdiscounting_choice_brms model
#'
#' Quantile credible intervals on the report-space-transformed draws;
#' columns `term`, `estimate`, `conf.low`, `conf.high`, `level`. `beta0`
#' stays on the identity scale across report spaces, matching
#' `tidy.beezdiscounting_choice()`.
#'
#' @param object A `beezdiscounting_choice_brms` object.
#' @param parm Optional terms: display names (`"k:(Intercept)"`, `"gamma"`,
#'   `"beta0"`) or TMB coefficient names (`"beta_k"`, `"log_gamma"`).
#' @param level Credible level (default 0.95).
#' @param report_space Reporting scale.
#' @param ... Unused.
#' @return A tibble.
#' @export
confint.beezdiscounting_choice_brms <- function(
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

  specs <- list(
    list(var = "b_logk_Intercept", term = "k:(Intercept)",
      tmb = "beta_k", transform = TRUE),
    list(var = "b_loggamma_Intercept", term = "gamma",
      tmb = "log_gamma", transform = TRUE)
  )
  if (isTRUE(object$param_info$intercept)) {
    specs <- c(specs, list(
      list(var = "b_b0_Intercept", term = "beta0",
        tmb = "beta0", transform = FALSE)
    ))
  }

  rows <- lapply(specs, function(s) {
    dj <- as.numeric(draws[, s$var])
    tdraws <- if (s$transform) .dd_brms_transform_draws(dj, to = to) else dj
    tibble::tibble(
      term = s$term,
      estimate = stats::median(tdraws),
      conf.low = unname(stats::quantile(tdraws, alpha2)),
      conf.high = unname(stats::quantile(tdraws, 1 - alpha2)),
      level = level
    )
  })
  out <- dplyr::bind_rows(rows)

  if (!is.null(parm)) {
    tmb_names <- vapply(specs, function(s) s$tmb, character(1))
    keep <- out$term %in% parm | tmb_names %in% parm
    out <- out[keep, , drop = FALSE]
  }
  out
}

#' @export
coef.beezdiscounting_choice_brms <- function(object, ...) {
  object$model$coefficients
}

#' @export
nobs.beezdiscounting_choice_brms <- function(object, ...) {
  object$param_info$n_obs
}

#' @export
logLik.beezdiscounting_choice_brms <- function(object, ...) {
  stop(
    "logLik() is not defined for Bayesian fits; use the stored $loo.",
    call. = FALSE
  )
}

#' Tidy a beezdiscounting_choice_brms model
#'
#' The 8-column dd contract; `gamma` (choice sensitivity) and `k` rows are
#' report-space transformed per draw; the optional `beta0` bias stays on
#' the logit scale.
#'
#' @param x A `beezdiscounting_choice_brms` object.
#' @param effects `"fixed"`, `"ran_pars"`, or both.
#' @param report_space Reporting scale.
#' @param ... Unused.
#' @return A tibble.
#' @export
tidy.beezdiscounting_choice_brms <- function(
  x,
  effects = c("fixed", "ran_pars"),
  report_space = c("natural", "log10", "internal", "log"),
  ...
) {
  effects <- match.arg(effects, several.ok = TRUE)
  report_space <- match.arg(report_space)

  result <- tibble::tibble()
  if ("fixed" %in% effects) {
    result <- dplyr::bind_rows(result, .dd_brms_choice_coef_table(x, report_space))
  }
  if ("ran_pars" %in% effects) {
    vc <- x$model$variance_components
    result <- dplyr::bind_rows(result, tibble::tibble(
      term = vc$Component,
      estimate = vc$Estimate,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      component = "variance",
      estimate_scale = vc$Scale,
      term_display = vc$Component
    ))
  }
  result
}

#' Glance at a beezdiscounting_choice_brms model
#'
#' @param x A `beezdiscounting_choice_brms` object.
#' @param ... Unused.
#' @return A one-row tibble.
#' @export
glance.beezdiscounting_choice_brms <- function(x, ...) {
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
    model_class = "beezdiscounting_choice_brms",
    backend = "brms",
    mode = "structural",
    equation = x$param_info$equation,
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

#' Predict P(LL) from a beezdiscounting_choice_brms model
#'
#' @param object A `beezdiscounting_choice_brms` object.
#' @param newdata Optional trial-level data (needs `rel` and `delay`, or
#'   `ss`/`ll` columns from which `rel` is derived).
#' @param type `"response"` (P(choosing the larger-later)) or
#'   `"parameters"` (per-subject k summaries).
#' @param level `"subject"` or `"population"`.
#' @param probs Interval probabilities.
#' @param draws Return the draws matrix instead.
#' @param ... Unused.
#' @return A tibble (or draws matrix).
#' @export
predict.beezdiscounting_choice_brms <- function(
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
  if (!is.null(newdata) && !"rel" %in% names(newdata)) {
    if (all(c("ss_amount", "ll_amount") %in% names(newdata))) {
      newdata$rel <- newdata$ll_amount / newdata$ss_amount
    } else {
      stop(
        "`newdata` needs a `rel` column (or `ss_amount` and `ll_amount`).",
        call. = FALSE
      )
    }
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
print.beezdiscounting_choice_brms <- function(x, digits = 4, ...) {
  cat("Bayesian Structural Choice Discounting Model (brms)\n")
  cat(strrep("=", 50), "\n")
  cat("Equation:", x$param_info$equation, "\n")
  cat(
    "Subjects:", x$param_info$n_subjects,
    " Trials:", x$param_info$n_obs, "\n"
  )
  cat("Converged:", ifelse(isTRUE(x$converged), "Yes", "No"), "\n")
  cat("\nCoefficients (posterior medians, natural scale):\n")
  ct <- .dd_brms_choice_coef_table(x, "natural")
  print(
    as.data.frame(ct[, c("term_display", "estimate", "std.error")]),
    digits = digits, row.names = FALSE
  )
  invisible(x)
}

#' Summarize a beezdiscounting_choice_brms model
#'
#' @param object A `beezdiscounting_choice_brms` object.
#' @param report_space Reporting scale.
#' @param ... Unused.
#' @return A `summary.beezdiscounting_choice_brms` list.
#' @export
summary.beezdiscounting_choice_brms <- function(
  object,
  report_space = c("natural", "log10", "internal", "log"),
  ...
) {
  report_space <- match.arg(report_space)
  structure(
    list(
      equation = object$param_info$equation,
      backend = "brms",
      converged = object$converged,
      n_subjects = object$param_info$n_subjects,
      nobs = object$param_info$n_obs,
      coefficients = .dd_brms_choice_coef_table(object, report_space),
      variance_components = object$model$variance_components,
      mcmc_info = object$mcmc_info,
      report_space = report_space
    ),
    class = "summary.beezdiscounting_choice_brms"
  )
}

#' @export
print.summary.beezdiscounting_choice_brms <- function(x, digits = 4, ...) {
  cat("\nBayesian Structural Choice Discounting Model Summary (brms)\n")
  cat(strrep("=", 50), "\n\n")
  cat("Equation:", x$equation, " Backend:", x$backend, "\n")
  cat("Converged:", ifelse(isTRUE(x$converged), "Yes", "No"), "\n")
  cat("Subjects:", x$n_subjects, " Trials:", x$nobs, "\n\n")
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
