#' Methods for linearized Mazur fits
#'
#' S3 methods for objects of class `beezdiscounting_linear` returned by
#' [fit_dd_linear()].
#' All quantities are closed form; nothing is re-optimized.
#'
#' `confint()` uses `parm` to SELECT WHICH INTERVAL is returned, not to filter
#' parameter names as in [stats::confint()].
#' `logLik()` reports the raw-scale (Jacobian-corrected) log-likelihood by
#' default, which makes it comparable
#' with this package's Gaussian NLS and TMB log-likelihoods.
#'
#' @param object,x A `beezdiscounting_linear` fit (for `print.summary()`, the
#'   object returned by `summary()`).
#' @param parm Which interval `confint()` returns: `"population"` (default) is
#'   an ANOVA-style t-interval on each condition mean, using the between-unit
#'   mean square with `N - C` degrees of freedom; `"subject"` gives the
#'   per-subject t-intervals on ln k. It does not select parameter names.
#' @param level For `confint()`, the confidence level (default `0.95`; note
#'   that `confint(parm = "subject")` uses this, not the fit's `conf_level`).
#'   For `logLik()`, which likelihood to return: `"population"` (the
#'   random-effects MLE, `df` = C + 2) or `"subject"` (the sum of the per-unit
#'   log-likelihoods, `df` = 2 per unit).
#' @param scale For `logLik()`: `"raw"` (default) is Jacobian-corrected back to
#'   the indifference-point scale; `"transformed"` is on the linearized scale.
#' @param effects For `tidy()`: `"subject"` (default) returns the per-unit
#'   table; `"population"` returns the random-effects MLEs `mu_<level>`,
#'   `sigma2` and `g`.
#' @param ... Unused; present for S3 generic consistency.
#' @return `print()` and `print.summary()` return their input invisibly;
#'   `summary()` returns an object of class `summary.beezdiscounting_linear`;
#'   `tidy()`, `glance()` and `augment()` return tibbles; `coef()` returns the
#'   named vector of condition means `mu`; `confint()` returns a two-column
#'   matrix of lower/upper bounds; `nobs()` returns the number of usable
#'   transformed observations; `logLik()` returns a `"logLik"` object with `df`
#'   and `nobs` attributes. Methods needing the random-effects component error
#'   when the design was unbalanced and `re` is `NULL`.
#' @seealso [fit_dd_linear()],
#'   [anova.beezdiscounting_linear()]
#' @name beezdiscounting_linear-methods
NULL

#' @rdname beezdiscounting_linear-methods
#' @export
print.beezdiscounting_linear <- function(x, ...) {
  cat("Linearized Mazur discounting fit (Hinds et al., 2026)\n")
  cat(sprintf(
    "  %d units, %d delays each%s\n",
    nrow(x$subjects),
    if (is.null(x$re)) NA_integer_ else x$re$n_delays,
    if (is.null(x$design$factor)) "" else sprintf(", factor: %s", x$design$factor)
  ))
  cat(sprintf(
    "  boundary = \"%s\" (eps = %g): %d point(s) at 0/1, %d clamped, %d dropped\n",
    x$transform$boundary, x$transform$eps, x$transform$n_boundary,
    x$transform$n_clamped, x$transform$n_dropped
  ))
  if (!is.null(x$re)) {
    cat("  Population ln k (mu):\n")
    print(round(x$re$mu, 3))
    cat(sprintf(
      "  sigma2 = %.3f, g = %.3f%s, logLik(raw) = %.2f\n",
      x$re$sigma2, x$re$g,
      if (x$re$g_zero) " [g-hat = 0 branch]" else "",
      x$re$loglik_raw
    ))
  } else {
    cat("  Random-effects component: not fitted (unbalanced design)\n")
  }
  invisible(x)
}

#' @rdname beezdiscounting_linear-methods
#' @export
summary.beezdiscounting_linear <- function(object, ...) {
  structure(
    list(
      fit = object,
      subjects = object$subjects,
      population = if (is.null(object$re)) NULL else tidy(object, effects = "population")
    ),
    class = "summary.beezdiscounting_linear"
  )
}

#' @rdname beezdiscounting_linear-methods
#' @export
print.summary.beezdiscounting_linear <- function(x, ...) {
  print(x$fit)
  cat("\nPer-unit ln k (first 10):\n")
  cols <- c("id", "condition", "logk", "se", "ci_lo", "ci_hi", "n_boundary")
  print(utils::head(x$subjects[, cols], 10))
  if (!is.null(x$population)) {
    cat("\nPopulation parameters:\n")
    print(x$population)
  }
  invisible(x)
}

#' @rdname beezdiscounting_linear-methods
#' @importFrom generics tidy
#' @export
tidy.beezdiscounting_linear <- function(x, effects = c("subject", "population"), ...) {
  effects <- match.arg(effects)
  if (effects == "subject") {
    return(x$subjects)
  }
  if (is.null(x$re)) {
    cli::cli_abort("Random-effects component was not fitted (unbalanced design).")
  }
  tibble::tibble(
    term = c(paste0("mu_", names(x$re$mu)), "sigma2", "g"),
    estimate = c(unname(x$re$mu), x$re$sigma2, x$re$g)
  )
}

#' @rdname beezdiscounting_linear-methods
#' @importFrom generics glance
#' @export
glance.beezdiscounting_linear <- function(x, ...) {
  ll <- if (is.null(x$re)) NA_real_ else x$re$loglik_raw
  n_par <- if (is.null(x$re)) NA_integer_ else x$re$n_par
  n <- nobs(x)
  tibble::tibble(
    model_class = "beezdiscounting_linear",
    backend = "closed_form",
    equation = "mazur_linearized",
    nobs = n,
    n_units = nrow(x$subjects),
    n_conditions = length(x$design$levels),
    n_par = n_par,
    sigma2 = if (is.null(x$re)) NA_real_ else x$re$sigma2,
    g = if (is.null(x$re)) NA_real_ else x$re$g,
    logLik = ll,
    AIC = -2 * ll + 2 * n_par,
    BIC = -2 * ll + log(n) * n_par
  )
}

#' @rdname beezdiscounting_linear-methods
#' @export
coef.beezdiscounting_linear <- function(object, ...) {
  if (is.null(object$re)) {
    cli::cli_abort("Random-effects component was not fitted.")
  }
  object$re$mu
}

#' @rdname beezdiscounting_linear-methods
#' @export
confint.beezdiscounting_linear <- function(object, parm = c("population", "subject"),
                                           level = 0.95, ...) {
  parm <- match.arg(parm)
  if (parm == "subject") {
    tq <- stats::qt(1 - (1 - level) / 2, df = object$subjects$df)
    out <- cbind(
      object$subjects$logk - tq * object$subjects$se,
      object$subjects$logk + tq * object$subjects$se
    )
    rownames(out) <- as.character(object$subjects$id)
  } else {
    re <- object$re
    if (is.null(re)) {
      cli::cli_abort("Random-effects component was not fitted.")
    }
    if (re$n_units - length(re$mu) < 1L) {
      cli::cli_abort(paste0(
        "population intervals need at least one condition with two or more units ",
        "(N - C = 0)"
      ))
    }
    # ANOVA-style interval on the per-subject ln k means (the paper gives no population
    # interval). Between-subject mean square with N - C df: SSR_{Z|X} / (T (N - C)).
    # The ML divisor N (= sigma2 (g+1)/T) would understate the SE by sqrt((N - C)/N).
    ms_between <- re$ssr_zx / (re$n_delays * (re$n_units - length(re$mu)))
    se <- sqrt(ms_between / as.numeric(re$n_per_condition[names(re$mu)]))
    tq <- stats::qt(1 - (1 - level) / 2, df = re$n_units - length(re$mu))
    out <- cbind(re$mu - tq * se, re$mu + tq * se)
    rownames(out) <- names(re$mu)
  }
  colnames(out) <- paste0(100 * c((1 - level) / 2, 1 - (1 - level) / 2), " %")
  out
}

#' @rdname beezdiscounting_linear-methods
#' @importFrom generics augment
#' @export
augment.beezdiscounting_linear <- function(x, ...) {
  d <- x$data
  k <- x$subjects$k[match(as.character(d$id), as.character(x$subjects$id))]
  d$.fitted <- 1 / (1 + k * d$x)
  d$.resid <- d$y - d$.fitted
  tibble::as_tibble(d)
}

#' @rdname beezdiscounting_linear-methods
#' @export
nobs.beezdiscounting_linear <- function(object, ...) {
  sum(is.finite(object$data$y_lin))
}

#' @rdname beezdiscounting_linear-methods
#' @export
logLik.beezdiscounting_linear <- function(object, level = c("population", "subject"),
                                          scale = c("raw", "transformed"), ...) {
  level <- match.arg(level)
  scale <- match.arg(scale)
  if (level == "population") {
    if (is.null(object$re)) {
      cli::cli_abort("Random-effects component was not fitted.")
    }
    val <- if (scale == "raw") object$re$loglik_raw else object$re$loglik_y
    df <- object$re$n_par
  } else {
    val <- sum(if (scale == "raw") object$subjects$loglik_raw else object$subjects$loglik_y)
    df <- 2L * nrow(object$subjects)
  }
  structure(val, df = df, nobs = nobs(object), class = "logLik")
}

#' F-test for condition means in a linearized Mazur fit (Hinds et al. 2026, Prop. 3.5)
#'
#' The paper's F-test is exact when the random-effects variance estimate
#' `g-hat > 0` (Prop. 3.5's assumption); when `g-hat = 0` the statistic is still
#' reported, with a warning.
#'
#' @param object A `beezdiscounting_linear` fit with a factor.
#' @param hypothesis `NULL` (all condition means equal) or a list of character
#'   vectors, each naming levels constrained equal under H0.
#' @param pairwise If `TRUE`, test every pair of levels (uncorrected p-values).
#' @param ... Unused; present for S3 generic consistency.
#' @return Tibble of class `beezdiscounting_linear_anova`.
#' @seealso [beezdiscounting_linear-methods] for the other S3 methods.
#' @export
anova.beezdiscounting_linear <- function(object, hypothesis = NULL, pairwise = FALSE, ...) {
  re <- object$re
  if (is.null(re)) {
    cli::cli_abort("Random-effects component was not fitted; {.code anova()} is unavailable.")
  }
  if (length(object$design$levels) < 2L) {
    cli::cli_abort("No conditions to compare: the fit has a single condition.")
  }
  if (pairwise) {
    pairs <- utils::combn(object$design$levels, 2, simplify = FALSE)
    out <- do.call(rbind, lapply(pairs, function(p) .dd_lin_ftest(re, list(p))))
  } else {
    out <- .dd_lin_ftest(re, hypothesis)
  }
  class(out) <- c("beezdiscounting_linear_anova", class(out))
  out
}
