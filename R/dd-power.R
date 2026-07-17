# ==============================================================================
# Monte Carlo power analysis for between-subject discounting designs
#
# power_discounting() simulates two-condition between-subject indifference-
# point data with simulate_dd_ip(), refits each replicate with fit_dd_tmb(),
# and estimates power as the proportion of usable fits whose Wald test on the
# log-k condition contrast rejects at the nominal alpha. find_n_discounting()
# is a bisection-search wrapper. Mirrors beezdemand::power_demand(); the
# within- vs between-subject asymmetry is inherited from each package's
# simulator by design.
# ==============================================================================

# Simulator-aligned design defaults (mirrors simulate_dd_ip()).
.dd_power_design_defaults <- function() {
  list(
    delays = c(7, 30, 180, 365, 730, 1460, 2920),
    log_k_pop = log(0.01),
    sigma_u = 0.6,
    phi = 10,
    sigma_e = 0.1
  )
}

#' Wilson score interval for a binomial proportion
#' @keywords internal
#' @noRd
.dd_power_wilson_ci <- function(x, n, level = 0.95) {
  if (n == 0) {
    return(c(NA_real_, NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- (z / denom) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, center - half), min(1, center + half))
}

#' Classify one Monte Carlo replicate and compute its Wald verdicts
#'
#' Pure function so the exclusion rules are unit-testable: a replicate enters
#' the power denominator only when the fit converged, the Hessian is positive
#' definite, and the target-term SE is finite. Anything else is surfaced via
#' `status` and its hit indicators are `NA` -- never counted as "no effect
#' detected".
#'
#' The Wald statistic is referred to a t distribution with `df` degrees of
#' freedom rather than the asymptotic normal: with plug-in variance
#' estimates, the null statistics are t-like at study-relevant N (the
#' asymptotic z-test was empirically anticonservative in the sibling
#' beezdemand calibration battery). `df = Inf` recovers the z-test.
#' @keywords internal
#' @noRd
.dd_power_rep_row <- function(
  converged,
  hessian_pd,
  estimate,
  se,
  alpha,
  df = Inf
) {
  status <- if (!isTRUE(converged)) {
    "nonconverged"
  } else if (!isTRUE(hessian_pd)) {
    "hessian_not_pd"
  } else if (
    !is.numeric(estimate) ||
      !is.finite(estimate) ||
      !is.numeric(se) ||
      !is.finite(se) ||
      se <= 0
  ) {
    "se_unusable"
  } else {
    "ok"
  }

  est_out <- if (is.numeric(estimate) && length(estimate) == 1) {
    as.numeric(estimate)
  } else {
    NA_real_
  }
  se_out <- if (is.numeric(se) && length(se) == 1) as.numeric(se) else NA_real_

  if (status != "ok") {
    return(list(
      status = status,
      estimate = est_out,
      se = se_out,
      statistic = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      hit_p = NA,
      hit_ci = NA
    ))
  }

  z <- estimate / se
  p <- 2 * stats::pt(-abs(z), df = df)
  crit <- stats::qt(1 - alpha / 2, df = df)
  ci_lower <- estimate - crit * se
  ci_upper <- estimate + crit * se
  list(
    status = status,
    estimate = est_out,
    se = se_out,
    statistic = z,
    p_value = p,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    hit_p = p < alpha,
    hit_ci = ci_lower > 0 || ci_upper < 0
  )
}

#' Validate the single-delta effect specification
#' @keywords internal
#' @noRd
.dd_power_validate_effect <- function(effect) {
  if (!is.list(effect)) {
    cli::cli_abort("{.arg effect} must be a list.")
  }
  .dd_power_require_named_unique(effect, "effect")
  bad <- setdiff(names(effect), "delta_k")
  if (length(bad) > 0) {
    cli::cli_abort(
      "{.arg effect} may only contain {.val delta_k}, not {.val {bad}}."
    )
  }
  if (is.null(effect$delta_k)) {
    cli::cli_abort(
      "Supply exactly one effect: {.code effect = list(delta_k = <shift on log k>)}."
    )
  }
  delta <- effect$delta_k
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta)) {
    cli::cli_abort(
      "`effect$delta_k` must be a single finite number (0 is allowed)."
    )
  }
  list(name = "delta_k", delta = as.numeric(delta))
}

#' Validate and merge the design list over simulator defaults
#' @keywords internal
#' @noRd
.dd_power_validate_design <- function(design) {
  defaults <- .dd_power_design_defaults()
  if (!is.list(design)) {
    cli::cli_abort("{.arg design} must be a list.")
  }
  .dd_power_require_named_unique(design, "design")
  bad <- setdiff(names(design), names(defaults))
  if (length(bad) > 0) {
    cli::cli_abort(c(
      "Unknown {.arg design} elements: {.val {bad}}.",
      "i" = "Allowed elements: {.val {names(defaults)}}."
    ))
  }
  merged <- utils::modifyList(defaults, design)
  for (nm in setdiff(names(defaults), "delays")) {
    v <- merged[[nm]]
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v)) {
      cli::cli_abort("`design${nm}` must be a single finite number.")
    }
  }
  for (nm in c("sigma_u", "phi", "sigma_e")) {
    if (merged[[nm]] <= 0) {
      cli::cli_abort("`design${nm}` must be > 0.")
    }
  }
  delays <- merged$delays
  if (
    !is.numeric(delays) ||
      length(delays) < 2 ||
      any(!is.finite(delays)) ||
      any(delays < 0)
  ) {
    cli::cli_abort(
      "`design$delays` must be a finite non-negative numeric vector of
       length >= 2."
    )
  }
  merged
}

#' Require every element of a user-supplied list to be uniquely named
#' @keywords internal
#' @noRd
.dd_power_require_named_unique <- function(x, arg) {
  if (length(x) == 0) {
    return(invisible(NULL))
  }
  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) {
    cli::cli_abort(
      "Every element of {.arg {arg}} must be named (unnamed elements would be
       silently ignored)."
    )
  }
  if (anyDuplicated(nms)) {
    cli::cli_abort(
      "{.arg {arg}} contains duplicated names: {.val {unique(nms[duplicated(nms)])}}."
    )
  }
  invisible(NULL)
}

#' Validate the seed argument
#' @keywords internal
#' @noRd
.dd_power_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  if (
    !is.numeric(seed) ||
      length(seed) != 1 ||
      !is.finite(seed) ||
      seed != round(seed) ||
      abs(seed) > .Machine$integer.max
  ) {
    cli::cli_abort(
      "{.arg seed} must be NULL or a single whole number representable as an
       integer."
    )
  }
  invisible(NULL)
}

#' Validate the df argument, resolving NULL to the design-based default
#' @keywords internal
#' @noRd
.dd_power_validate_df <- function(df, default_df) {
  if (is.null(df)) {
    if (!is.null(default_df) && default_df <= 0) {
      cli::cli_abort(
        "{.arg n_subjects} is too small for the design-based default {.arg df}
         ({default_df}); increase {.arg n_subjects} or supply {.arg df}."
      )
    }
    return(default_df)
  }
  if (!is.numeric(df) || length(df) != 1 || is.na(df) || df <= 0) {
    cli::cli_abort(
      "{.arg df} must be NULL (design-based default), a positive number, or Inf."
    )
  }
  df
}

#' Shared scalar-argument validation for the power engine
#' @keywords internal
#' @noRd
.dd_power_validate_scalars <- function(n_subjects, n_sim, alpha) {
  if (
    !is.numeric(n_subjects) ||
      length(n_subjects) != 1 ||
      !is.finite(n_subjects) ||
      n_subjects != round(n_subjects) ||
      n_subjects < 2
  ) {
    cli::cli_abort("{.arg n_subjects} must be a single whole number >= 2.")
  }
  if (
    !is.numeric(n_sim) ||
      length(n_sim) != 1 ||
      !is.finite(n_sim) ||
      n_sim != round(n_sim) ||
      n_sim < 1
  ) {
    cli::cli_abort("{.arg n_sim} must be a single whole number >= 1.")
  }
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1 ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    cli::cli_abort(
      "{.arg alpha} must be a single number strictly between 0 and 1."
    )
  }
  invisible(NULL)
}

#' Reject fit arguments that the power engine sets itself
#' @keywords internal
#' @noRd
.dd_power_check_fit_args <- function(fit_args) {
  if (
    length(fit_args) > 0 &&
      (is.null(names(fit_args)) || any(!nzchar(names(fit_args))))
  ) {
    cli::cli_abort(
      "All arguments passed via {.arg ...} must be named (they are forwarded
       to the fitting function)."
    )
  }
  reserved <- c("data", "y_var", "x_var", "id_var", "factors")
  bad <- intersect(names(fit_args), reserved)
  if (length(bad) > 0) {
    cli::cli_abort(
      "Arguments {.val {bad}} are set by the power engine and cannot be
       overridden."
    )
  }
  fit_args
}

#' Run the discounting Monte Carlo replicate loop
#'
#' Each replicate simulates a two-condition between-subject dataset with
#' `simulate_dd_ip()`, refits it with [fit_dd_tmb()], and extracts the log-k
#' condition contrast on the estimation (natural log) scale via
#' `tidy(fit, report_space = "internal")`. Replicate-level errors are caught
#' and recorded as `status = "error"`, never propagated.
#' @keywords internal
#' @noRd
.dd_power_replicates <- function(
  n_subjects,
  delta,
  design,
  n_sim,
  alpha,
  df,
  equation,
  family,
  random_effects,
  multi_start,
  fit_args,
  sim_offset = 0L,
  verbose = FALSE
) {
  target_term <- "k:conditionC2"
  rows <- vector("list", n_sim)
  # Force treatment coding for the duration of the loop so the target term
  # name and estimand (C2 - C1 on the log-k scale) survive user-level
  # options(contrasts = ...) such as contr.sum. The factor-attribute route
  # is not robust: internal factor handling in the fitter may rebuild the
  # factor and drop a per-factor contrasts attribute.
  old_opts <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old_opts), add = TRUE)
  if (verbose) {
    cli::cli_progress_bar(
      "Monte Carlo replicates (n_subjects = {n_subjects})",
      total = n_sim
    )
  }

  for (i in seq_len(n_sim)) {
    rows[[i]] <- tryCatch(
      {
        sim <- simulate_dd_ip(
          n_subjects = n_subjects,
          delays = design$delays,
          log_k_pop = design$log_k_pop,
          sigma_u = design$sigma_u,
          phi = design$phi,
          sigma_e = design$sigma_e,
          family = family,
          equation = equation,
          n_conditions = 2,
          delta_k = c(0, delta),
          seed = NULL
        )
        fit <- suppressWarnings(suppressMessages(do.call(
          fit_dd_tmb,
          c(
            list(
              data = sim,
              y_var = "y",
              x_var = "x",
              id_var = "id",
              equation = equation,
              family = family,
              factors = "condition",
              random_effects = random_effects,
              multi_start = multi_start,
              verbose = 0
            ),
            fit_args
          )
        )))

        est <- NA_real_
        se <- NA_real_
        extract_msg <- NA_character_
        if (!isTRUE(fit$converged)) {
          opt_msg <- tryCatch(
            as.character(fit$opt$message)[1],
            error = function(e) NA_character_
          )
          if (length(opt_msg) == 1 && !is.na(opt_msg) && nzchar(opt_msg)) {
            extract_msg <- opt_msg
          }
        }
        if (isTRUE(fit$converged)) {
          td <- tryCatch(
            tidy(fit, effects = "fixed", report_space = "internal"),
            error = function(e) {
              extract_msg <<- paste0("tidy() failed: ", conditionMessage(e))
              NULL
            }
          )
          if (!is.null(td)) {
            hit_row <- td[td$term == target_term, , drop = FALSE]
            if (nrow(hit_row) == 1) {
              est <- hit_row$estimate
              se <- hit_row$std.error
            } else {
              extract_msg <- sprintf(
                "target term '%s' matched %d rows in tidy() output",
                target_term,
                nrow(hit_row)
              )
            }
          }
        }

        rep_out <- .dd_power_rep_row(
          converged = isTRUE(fit$converged),
          hessian_pd = isTRUE(fit$hessian_pd),
          estimate = est,
          se = se,
          alpha = alpha,
          df = df
        )
        c(
          list(
            sim = i + sim_offset,
            converged = isTRUE(fit$converged),
            hessian_pd = isTRUE(fit$hessian_pd)
          ),
          rep_out,
          list(message = extract_msg)
        )
      },
      error = function(e) {
        row <- c(
          list(sim = i + sim_offset, converged = FALSE, hessian_pd = FALSE),
          .dd_power_rep_row(FALSE, FALSE, NA_real_, NA_real_, alpha, df),
          list(message = conditionMessage(e))
        )
        # A caught execution error is distinct from ordinary optimizer
        # nonconvergence (whose optimizer message is surfaced above).
        row$status <- "error"
        row
      }
    )
    if (verbose) cli::cli_progress_update()
  }
  if (verbose) {
    cli::cli_progress_done()
  }

  out <- dplyr::bind_rows(rows)
  out[, c(
    "sim",
    "status",
    "converged",
    "hessian_pd",
    "estimate",
    "se",
    "statistic",
    "p_value",
    "ci_lower",
    "ci_upper",
    "hit_p",
    "hit_ci",
    "message"
  )]
}

#' Summarize a replicate table into power estimates and diagnostics
#' @keywords internal
#' @noRd
.dd_power_summarize <- function(replicates, mc_ci_level = 0.95) {
  ok <- replicates$status == "ok"
  n_used <- sum(ok)
  hits_ci <- sum(replicates$hit_ci[ok])
  hits_p <- sum(replicates$hit_p[ok])
  list(
    power = if (n_used > 0) hits_ci / n_used else NA_real_,
    power_mc_ci = .dd_power_wilson_ci(hits_ci, n_used, mc_ci_level),
    hit_rate_p = if (n_used > 0) hits_p / n_used else NA_real_,
    hit_rate_ci = if (n_used > 0) hits_ci / n_used else NA_real_,
    n_sim = nrow(replicates),
    n_converged = sum(replicates$converged),
    n_hessian_pd = sum(replicates$hessian_pd),
    n_used = n_used
  )
}

#' Warn when too few replicates produced usable fits
#' @keywords internal
#' @noRd
.dd_power_usable_fraction_warn <- function(n_used, n_sim, threshold = 0.95) {
  if (n_sim > 0 && n_used / n_sim < threshold) {
    cli::cli_warn(c(
      "!" = "Only {n_used}/{n_sim} replicate{?s} produced usable fits
             (converged, positive-definite Hessian, finite SE).",
      "i" = "The reported power is conditional on a usable fit and can be
             selected when convergence depends on the realized data. Inspect
             {.code $replicates$status} and consider a larger {.arg n_subjects}
             or {.code multi_start = TRUE}."
    ))
  }
  invisible(NULL)
}

#' Monte Carlo power analysis for between-subject discounting designs
#'
#' @description
#' Estimates statistical power to detect a between-subject difference in the
#' discount rate (a shift `delta_k` on natural-log k) by simulation: each
#' replicate (1) simulates a two-condition between-subject indifference-point
#' dataset with [simulate_dd_ip()] under assumed population parameters plus
#' the effect, (2) refits it with [fit_dd_tmb()], and (3) tests the log-k
#' condition contrast with a Wald test at level `alpha`, referred to a t
#' distribution with `df` degrees of freedom (see the `df` argument). Power
#' is the proportion of *usable* fits (converged, positive-definite Hessian,
#' finite standard error) that reject.
#'
#' Because the power estimate is a proportion from finitely many replicates,
#' it is reported with a Wilson score confidence interval (`power_mc_ci`).
#' Both a p-value verdict (`p < alpha`) and a confidence-interval verdict
#' (Wald CI excludes 0) are recorded per replicate; they use the same
#' standard error and reference distribution, so they coincide by
#' construction, and both rates are returned.
#'
#' This mirrors `beezdemand::power_demand()`. The design asymmetry is
#' intentional: this package's simulator models a *between-subject* condition
#' (subjects are split across conditions round-robin), while beezdemand's
#' models a within-subject condition.
#'
#' @param n_subjects Total number of simulated subjects per replicate, split
#'   across the two conditions round-robin (even numbers give equal groups).
#' @param effect Named list supplying `delta_k`: the true condition-2 shift
#'   on natural-log k. `0` is allowed (useful for Type I error checks).
#'   E.g. `delta_k = log(2)` means condition 2's k is twice condition 1's.
#' @param design Named list of data-generating settings, merged over the
#'   simulator defaults: `delays` (vector), `log_k_pop`, `sigma_u` (subject
#'   SD on log k), `phi` (SLT-beta precision; used when `family = "sltb"`),
#'   and `sigma_e` (residual SD; used when `family = "gaussian"`).
#' @param n_sim Number of Monte Carlo replicates. 500 (default) is suitable
#'   for interactive exploration; use 2000+ for grant-quality precision (see
#'   `vignette("power-analysis")`).
#' @param alpha Nominal two-sided test level.
#' @param df Degrees of freedom for the Wald test's t reference
#'   distribution. `NULL` (default) uses `n_subjects - 2`, the two-sample df
#'   of the between-subject design. This is an *empirically calibrated*
#'   small-sample correction, not a model-derived df (the TMB fit has no
#'   exact t sampling theory); it passes the package's Type I calibration
#'   battery, while the asymptotic z-test (`df = Inf`) is anticonservative
#'   at study-relevant sample sizes.
#' @param seed Optional integer seed; identical seeds give identical
#'   results. The caller's RNG state is restored on exit.
#' @param equation Discounting function used for BOTH simulation and
#'   refitting (they are always matched): `"mazur"` (default) or
#'   `"exponential"`. The two-parameter Green-Myerson and Rachlin forms are
#'   out of scope in this version.
#' @param family Observation family used for both simulation and refitting:
#'   `"sltb"` (default) or `"gaussian"`.
#' @param random_effects Random-effects formula passed to [fit_dd_tmb()].
#'   The default `k ~ 1` matches the simulator's data-generating process (a
#'   single subject random intercept on log k).
#' @param multi_start Passed to [fit_dd_tmb()]. Defaults to `FALSE` for
#'   speed; non-convergent replicates are excluded and surfaced rather than
#'   biasing the estimate.
#' @param verbose Logical; show a progress bar.
#' @param ... Additional arguments passed to [fit_dd_tmb()] (e.g.
#'   `tmb_control`).
#'
#' @return An object of class `beezdiscounting_power`: a list with
#'   \describe{
#'     \item{power}{Estimated power: proportion of usable replicates whose
#'       Wald CI excludes 0 (equal to `hit_rate_ci`). `NA` if no replicate
#'       was usable.}
#'     \item{power_mc_ci}{Wilson 95% confidence interval on `power`,
#'       reflecting Monte Carlo uncertainty from `n_used` replicates.}
#'     \item{hit_rate_p}{Proportion of usable replicates with `p < alpha`.}
#'     \item{hit_rate_ci}{Proportion of usable replicates whose Wald CI
#'       excludes 0 (the same decision rule as `hit_rate_p`, since both use
#'       the same SE and t reference; both reported).}
#'     \item{n_sim}{Total replicates attempted.}
#'     \item{n_converged}{Replicates whose fit converged.}
#'     \item{n_hessian_pd}{Replicates with a positive-definite Hessian.}
#'     \item{n_used}{Replicates entering the power denominator (converged,
#'       positive-definite Hessian, finite SE).}
#'     \item{alpha}{Nominal test level.}
#'     \item{df}{Degrees of freedom of the t reference distribution actually
#'       used (`n_subjects - 2` unless overridden).}
#'     \item{effect}{The validated effect specification (name and delta).}
#'     \item{target_term}{The tested coefficient (`"k:conditionC2"`).}
#'     \item{design}{The merged design list actually used.}
#'     \item{n_subjects}{As supplied.}
#'     \item{replicates}{Tibble with one row per replicate: `sim`, `status`
#'       (`"ok"`, `"nonconverged"`, `"hessian_not_pd"`, `"se_unusable"`,
#'       `"error"`), `converged`, `hessian_pd`, `estimate`, `se`,
#'       `statistic`, `p_value`, `ci_lower`, `ci_upper`, `hit_p`, `hit_ci`,
#'       and `message` (error text, if any). Estimates are on the
#'       natural-log-k scale of the simulated `delta_k`.}
#'     \item{seed}{As supplied.}
#'     \item{settings}{List of `equation`, `family`, `multi_start`, and the
#'       deparsed random-effects formula.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @details
#' A replicate whose fit fails (non-convergence, non-positive-definite
#' Hessian, unusable standard error, or an error) is excluded from the power
#' denominator and reported through the `n_*` counts and
#' `$replicates$status` -- it is never counted as "no effect detected", which
#' would bias power in an unpredictable direction. A warning is issued when
#' fewer than 95% of replicates are usable.
#'
#' The v1 scope is a single fixed-effect `delta_k` under the package's
#' existing between-subject simulator. Effects on `s` or `phi`, power for
#' derived measures (ED50, AUC), and arbitrary designs are out of scope; see
#' `vignette("power-analysis")`.
#'
#' @examples
#' \donttest{
#' # Quick exploratory run (use n_sim >= 500 for real planning)
#' res <- power_discounting(
#'   n_subjects = 40,
#'   effect = list(delta_k = log(2)),
#'   n_sim = 20, seed = 1, verbose = FALSE
#' )
#' print(res)
#' }
#'
#' @seealso [find_n_discounting()] to search for the smallest adequate
#'   sample size; [fit_dd_tmb()] for the model being refit;
#'   [simulate_dd_ip()] for the data-generating process.
#' @family power-analysis
#' @export
power_discounting <- function(
  n_subjects,
  effect = list(delta_k = NULL),
  design = list(),
  n_sim = 500,
  alpha = 0.05,
  df = NULL,
  seed = NULL,
  equation = c("mazur", "exponential"),
  family = c("sltb", "gaussian"),
  random_effects = k ~ 1,
  multi_start = FALSE,
  verbose = TRUE,
  ...
) {
  cl <- match.call()
  equation <- match.arg(equation)
  family <- match.arg(family)
  eff <- .dd_power_validate_effect(effect)
  design <- .dd_power_validate_design(design)
  .dd_power_validate_scalars(n_subjects, n_sim, alpha)
  .dd_power_validate_seed(seed)
  df <- .dd_power_validate_df(df, n_subjects - 2)
  fit_args <- .dd_power_check_fit_args(list(...))

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (
        exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      ) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
    set.seed(seed)
  }

  replicates <- .dd_power_replicates(
    n_subjects = n_subjects,
    delta = eff$delta,
    design = design,
    n_sim = n_sim,
    alpha = alpha,
    df = df,
    equation = equation,
    family = family,
    random_effects = random_effects,
    multi_start = multi_start,
    fit_args = fit_args,
    verbose = verbose
  )
  s <- .dd_power_summarize(replicates)
  .dd_power_usable_fraction_warn(s$n_used, s$n_sim)

  structure(
    c(
      s,
      list(
        alpha = alpha,
        df = df,
        effect = eff,
        target_term = "k:conditionC2",
        design = design,
        n_subjects = n_subjects,
        replicates = replicates,
        seed = seed,
        settings = list(
          equation = equation,
          family = family,
          multi_start = multi_start,
          random_effects = paste(deparse(random_effects), collapse = " ")
        ),
        call = cl
      )
    ),
    class = "beezdiscounting_power"
  )
}

#' @export
print.beezdiscounting_power <- function(x, ...) {
  cat("Monte Carlo power analysis (beezdiscounting)\n")
  cat(sprintf(
    "  Target: %s (%s = %.4g), two-sided alpha = %g, t reference (df = %g)\n",
    x$target_term,
    x$effect$name,
    x$effect$delta,
    x$alpha,
    x$df
  ))
  cat(sprintf(
    "  n_subjects = %d (between-subject, 2 conditions), n_sim = %d (converged %d, usable %d)\n",
    as.integer(x$n_subjects),
    as.integer(x$n_sim),
    as.integer(x$n_converged),
    as.integer(x$n_used)
  ))
  if (is.na(x$power)) {
    cat("  Power: NA (no usable fits)\n")
  } else {
    cat(sprintf(
      "  Power (CI-exclusion): %.3f [95%% MC CI %.3f, %.3f]\n",
      x$power,
      x$power_mc_ci[1],
      x$power_mc_ci[2]
    ))
    cat(sprintf("  p-value hit rate:     %.3f\n", x$hit_rate_p))
  }
  invisible(x)
}

# ==============================================================================
# Sample-size search
# ==============================================================================

#' Evaluate power at one N with adaptive replication at ambiguous results
#' @keywords internal
#' @noRd
.dd_power_eval_adaptive <- function(
  run_batch,
  n,
  target_power,
  n_sim,
  n_sim_max,
  verbose = FALSE
) {
  replicates <- run_batch(n, n_sim, 0L)
  repeat {
    s <- .dd_power_summarize(replicates)
    if (s$n_used == 0) {
      cli::cli_abort(c(
        "No usable fits at {.arg n_subjects} = {n} ({s$n_sim} replicates
         attempted).",
        "i" = "Power cannot be evaluated here; check the design or effect
               specification with {.fn power_discounting} directly."
      ))
    }
    ci <- s$power_mc_ci
    decision <- if (ci[1] >= target_power) {
      "above"
    } else if (ci[2] < target_power) {
      "below"
    } else {
      "ambiguous"
    }
    if (decision != "ambiguous" || s$n_sim >= n_sim_max) {
      break
    }
    batch <- min(n_sim, n_sim_max - s$n_sim)
    replicates <- dplyr::bind_rows(
      replicates,
      run_batch(n, batch, nrow(replicates))
    )
  }

  uncertain <- decision == "ambiguous"
  if (uncertain) {
    decision <- if (s$power >= target_power) {
      "ambiguous_above"
    } else {
      "ambiguous_below"
    }
  }
  if (verbose) {
    cli::cli_inform(
      "n_subjects = {n}: power {sprintf('%.3f', s$power)}
       [{sprintf('%.3f', ci[1])}, {sprintf('%.3f', ci[2])}]
       ({s$n_sim} sims) -> {decision}"
    )
  }
  list(
    row = tibble::tibble(
      n_subjects = n,
      n_sim_total = s$n_sim,
      n_used = s$n_used,
      usable_fraction = s$n_used / s$n_sim,
      power = s$power,
      ci_lower = ci[1],
      ci_upper = ci[2],
      decision = decision
    ),
    above = decision %in% c("above", "ambiguous_above"),
    uncertain = uncertain,
    power = s$power
  )
}

#' Shared bisection search over n_subjects
#' @keywords internal
#' @noRd
.dd_power_find_n_search <- function(
  run_batch,
  target_power,
  n_range,
  n_sim,
  n_sim_max,
  verbose
) {
  if (
    !is.numeric(n_range) ||
      length(n_range) != 2 ||
      any(!is.finite(n_range)) ||
      any(n_range != round(n_range)) ||
      n_range[1] < 2 ||
      n_range[1] >= n_range[2]
  ) {
    cli::cli_abort(
      "{.arg n_range} must be two whole numbers with 2 <= n_range[1] < n_range[2]."
    )
  }
  if (
    !is.numeric(target_power) ||
      length(target_power) != 1 ||
      !is.finite(target_power) ||
      target_power <= 0 ||
      target_power >= 1
  ) {
    cli::cli_abort(
      "{.arg target_power} must be a single number strictly between 0 and 1."
    )
  }

  lo <- as.integer(n_range[1])
  hi <- as.integer(n_range[2])
  evals <- list()
  any_uncertain <- FALSE

  eval_n <- function(n) {
    res <- .dd_power_eval_adaptive(
      run_batch,
      n,
      target_power,
      n_sim,
      n_sim_max,
      verbose = verbose
    )
    evals[[length(evals) + 1]] <<- res$row
    any_uncertain <<- any_uncertain || res$uncertain
    res
  }

  ev_hi <- eval_n(hi)
  if (!ev_hi$above) {
    cli::cli_abort(c(
      "Estimated power at the upper bound of {.arg n_range}
       ({.val {hi}} subjects) is {sprintf('%.3f', ev_hi$power)}, which does
       not reach the target of {target_power}.",
      "i" = "Increase {.arg n_range}, the effect size, or reconsider the design."
    ))
  }

  finish <- function(n, status, uncertain) {
    evaluations <- dplyr::bind_rows(evals)
    min_usable <- min(evaluations$usable_fraction)
    if (min_usable < 0.95) {
      cli::cli_warn(c(
        "!" = "At least one evaluated sample size had fewer than 95% usable
               fits (minimum usable fraction
               {sprintf('%.2f', min_usable)}).",
        "i" = "Power at those N is conditional on a usable fit; inspect
               {.code $evaluations} and consider {.code multi_start = TRUE}."
      ))
    }
    list(
      n = n,
      status = status,
      uncertain = uncertain,
      evaluations = evaluations
    )
  }

  ev_lo <- eval_n(lo)
  if (ev_lo$above) {
    return(finish(lo, "at_lower_bound", any_uncertain))
  }

  while (hi - lo > 1) {
    mid <- as.integer(floor((lo + hi) / 2))
    if (eval_n(mid)$above) hi <- mid else lo <- mid
  }

  # Confirmation pass: re-evaluate the selected N and its lower neighbor with
  # fresh replicates before claiming anything. Three outcomes:
  # - selected N re-confirms above and N - 1 below: "confirmed" (or
  #   "uncertain" if any decision along the way used a point estimate);
  # - selected N FAILS reconfirmation: the search evidence is contradicted;
  #   return n = NA with status "unresolved" rather than a number the run
  #   itself does not support;
  # - N - 1 also clears the target on reconfirmation: N reaches the target
  #   but may not be minimal -> "uncertain".
  conf_hi <- eval_n(hi)
  if (!conf_hi$above) {
    return(finish(NA_integer_, "unresolved", TRUE))
  }
  conf_lo_above <- if (hi - 1 >= n_range[1]) eval_n(hi - 1)$above else FALSE
  status <- if (conf_lo_above || any_uncertain) "uncertain" else "confirmed"

  finish(hi, status, any_uncertain || status != "confirmed")
}

#' Find the smallest sample size reaching a target power (discounting)
#'
#' @description
#' Bisection search over `n_subjects` for the smallest total N whose Monte
#' Carlo power estimate from [power_discounting()] reaches `target_power`.
#' The search accounts for Monte Carlo noise: at each evaluated N,
#' replicates are added in batches (up to `n_sim_max`) until the Wilson
#' interval for power lies wholly above or below the target; if it still
#' straddles the target at the cap, the decision falls back to the point
#' estimate and the result is flagged `uncertain`. The selected N and its
#' lower neighbor are then re-evaluated with fresh replicates before
#' minimality is claimed.
#'
#' The returned `n` is an *estimated minimum under Monte Carlo uncertainty*,
#' not an exact bound. For grant-quality reporting, rerun
#' [power_discounting()] at the returned `n` with a large `n_sim` (2000+)
#' and report that estimate with its Monte Carlo confidence interval.
#'
#' @param target_power Target power in (0, 1).
#' @inheritParams power_discounting
#' @param n_range Integer bracket `c(lower, upper)` of total N to search
#'   (`2 <= lower < upper`). The search errors (rather than extrapolating)
#'   if the target is not reached at `upper`.
#' @param n_sim Replicates per evaluation batch. Smaller than the
#'   [power_discounting()] default because several N values are evaluated;
#'   the adaptive rule adds batches where the verdict is close.
#' @param n_sim_max Maximum replicates per evaluated N (default `4 * n_sim`).
#' @param df Degrees of freedom for the Wald test's t reference. `NULL`
#'   (default) tracks the evaluated sample size as `n - 2`; a numeric value
#'   (or `Inf` for the asymptotic z-test) is used at every evaluated N.
#' @param verbose Logical; report each evaluation.
#' @param ... Additional arguments passed to [fit_dd_tmb()].
#'
#' @return An object of class `beezdiscounting_power_n`: a list with
#'   \describe{
#'     \item{n}{Estimated smallest total `n_subjects` reaching
#'       `target_power`; `NA` when the confirmation pass contradicted the
#'       search (`status = "unresolved"`).}
#'     \item{target_power}{As supplied.}
#'     \item{status}{`"confirmed"` (selected N re-confirmed above target and
#'       N - 1 below), `"uncertain"` (a decision relied on a point estimate,
#'       or N - 1 also cleared the target on reconfirmation so the returned
#'       N may not be minimal), `"unresolved"` (the selected N failed its
#'       own reconfirmation; `n` is `NA`), or `"at_lower_bound"` (the target
#'       was already met at `n_range[1]`; smaller N was not explored). These
#'       labels describe a heuristic Monte Carlo decision rule -- repeated
#'       looks at ordinary Wilson intervals across several N -- not a formal
#'       sequential error guarantee.}
#'     \item{uncertain}{Logical; `TRUE` when any search decision was made on
#'       a point estimate rather than a conclusive Wilson interval, or the
#'       status is not `"confirmed"`/`"at_lower_bound"`.}
#'     \item{evaluations}{Tibble of every evaluation: `n_subjects`,
#'       `n_sim_total`, `n_used`, `usable_fraction`, `power`, `ci_lower`,
#'       `ci_upper`, `decision`. A warning fires when any evaluation had
#'       fewer than 95% usable fits.}
#'     \item{alpha, df, effect, design, n_range, n_sim, n_sim_max, seed,
#'       settings, call}{Echoed inputs and effective settings.}
#'   }
#'
#' @examples
#' \donttest{
#' # Small search for demonstration (use larger n_sim for real planning)
#' res <- find_n_discounting(
#'   target_power = 0.8,
#'   effect = list(delta_k = log(4)),
#'   n_range = c(4, 40), n_sim = 30, seed = 1, verbose = FALSE
#' )
#' print(res)
#' }
#'
#' @seealso [power_discounting()] for the Monte Carlo engine.
#' @family power-analysis
#' @export
find_n_discounting <- function(
  target_power = 0.8,
  effect = list(delta_k = NULL),
  design = list(),
  n_range = c(6, 200),
  n_sim = 200,
  n_sim_max = 4 * n_sim,
  alpha = 0.05,
  df = NULL,
  seed = NULL,
  equation = c("mazur", "exponential"),
  family = c("sltb", "gaussian"),
  random_effects = k ~ 1,
  multi_start = FALSE,
  verbose = TRUE,
  ...
) {
  cl <- match.call()
  equation <- match.arg(equation)
  family <- match.arg(family)
  eff <- .dd_power_validate_effect(effect)
  design <- .dd_power_validate_design(design)
  .dd_power_validate_scalars(n_subjects = 2, n_sim = n_sim, alpha = alpha)
  .dd_power_validate_seed(seed)
  if (!is.null(df)) {
    .dd_power_validate_df(df, NULL)
  }
  if (
    !is.numeric(n_sim_max) ||
      length(n_sim_max) != 1 ||
      !is.finite(n_sim_max) ||
      n_sim_max != round(n_sim_max) ||
      n_sim_max < n_sim
  ) {
    cli::cli_abort(
      "{.arg n_sim_max} must be a single whole number >= {.arg n_sim}."
    )
  }
  fit_args <- .dd_power_check_fit_args(list(...))

  if (
    is.null(df) &&
      is.numeric(n_range) &&
      length(n_range) == 2 &&
      is.finite(n_range[1]) &&
      n_range[1] < 3
  ) {
    cli::cli_abort(
      "{.code n_range[1]} must be >= 3 when {.arg df} is NULL (the default
       df is n - 2)."
    )
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (
        exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      ) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
    set.seed(seed)
  }

  run_batch <- function(n, batch_size, sim_offset) {
    .dd_power_replicates(
      n_subjects = n,
      delta = eff$delta,
      design = design,
      n_sim = batch_size,
      alpha = alpha,
      df = if (is.null(df)) n - 2 else df,
      equation = equation,
      family = family,
      random_effects = random_effects,
      multi_start = multi_start,
      fit_args = fit_args,
      sim_offset = sim_offset,
      verbose = FALSE
    )
  }

  search <- .dd_power_find_n_search(
    run_batch,
    target_power,
    n_range,
    n_sim,
    n_sim_max,
    verbose
  )

  structure(
    c(
      search,
      list(
        target_power = target_power,
        alpha = alpha,
        df = df,
        effect = eff,
        design = design,
        n_range = n_range,
        n_sim = n_sim,
        n_sim_max = n_sim_max,
        seed = seed,
        settings = list(
          equation = equation,
          family = family,
          multi_start = multi_start,
          random_effects = paste(deparse(random_effects), collapse = " ")
        ),
        call = cl
      )
    ),
    class = "beezdiscounting_power_n"
  )
}

#' @export
print.beezdiscounting_power_n <- function(x, ...) {
  cat("Sample-size search (Monte Carlo power)\n")
  cat(sprintf(
    "  Target power %.2f for %s = %.4g at alpha = %g\n",
    x$target_power,
    x$effect$name,
    x$effect$delta,
    x$alpha
  ))
  if (is.na(x$n)) {
    cat("  No sample size confirmed (status: unresolved) --\n")
    cat("  the confirmation pass contradicted the search; rerun with a\n")
    cat("  larger n_sim / n_sim_max or inspect $evaluations.\n")
  } else {
    cat(sprintf(
      "  Estimated minimum n_subjects = %d (status: %s)\n",
      as.integer(x$n),
      x$status
    ))
    cat("  This is an estimated minimum under Monte Carlo uncertainty;\n")
    cat(
      "  rerun the power function at this N with a large n_sim to report it.\n"
    )
  }
  cat("\n  Evaluations:\n")
  print(as.data.frame(x$evaluations), row.names = FALSE)
  invisible(x)
}
