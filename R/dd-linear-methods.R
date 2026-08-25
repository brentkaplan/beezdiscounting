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

#' @export
coef.beezdiscounting_linear <- function(object, ...) {
  if (is.null(object$re)) {
    cli::cli_abort("Random-effects component was not fitted.")
  }
  object$re$mu
}

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

#' @importFrom generics augment
#' @export
augment.beezdiscounting_linear <- function(x, ...) {
  d <- x$data
  k <- x$subjects$k[match(as.character(d$id), as.character(x$subjects$id))]
  d$.fitted <- 1 / (1 + k * d$x)
  d$.resid <- d$y - d$.fitted
  tibble::as_tibble(d)
}

#' @export
nobs.beezdiscounting_linear <- function(object, ...) {
  sum(is.finite(object$data$y_lin))
}

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
#' @param object A `beezdiscounting_linear` fit with a factor.
#' @param hypothesis `NULL` (all condition means equal) or a list of character
#'   vectors, each naming levels constrained equal under H0.
#' @param pairwise If `TRUE`, test every pair of levels (uncorrected p-values).
#' @param ... Unused; present for S3 generic consistency.
#' @return Tibble of class `beezdiscounting_linear_anova`.
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
