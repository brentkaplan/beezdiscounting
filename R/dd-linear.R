#' Linearized-Mazur transform (Hinds et al. 2026, Eq. 6) with a boundary policy
#'
#' `y = ln(1/D - 1) - ln(t)` is undefined at `D = 0` and `D = 1`. The paper is
#' silent on those values; the package clamps them to `[eps, 1 - eps]` by
#' default (spec D1).
#' @param d numeric, normalized indifference points in `[0, 1]`.
#' @param t numeric, positive delays.
#' @param boundary one of `"clamp"`, `"drop"`, `"error"`.
#' @param eps clamp half-width.
#' @return list(y_lin, d_used, n_boundary, n_clamped, n_dropped, log_jac).
#' @keywords internal
#' @noRd
.dd_lin_transform <- function(d, t, boundary = c("clamp", "drop", "error"),
                              eps = 0.005) {
  boundary <- match.arg(boundary)
  if (any(!is.finite(t)) || any(t <= 0)) {
    cli::cli_abort("Delays {.arg t} must be finite and positive; ln(t) is undefined otherwise.")
  }
  if (any(!is.finite(d)) || any(d < 0 | d > 1)) {
    cli::cli_abort("Indifference points must lie in [0, 1] after scale conversion.")
  }
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0 || eps >= 0.5) {
    cli::cli_abort("{.arg eps} must be a single number in (0, 0.5).")
  }
  is_boundary <- d == 0 | d == 1
  n_boundary <- sum(is_boundary)
  d_used <- d
  n_clamped <- 0L
  n_dropped <- 0L
  if (boundary == "error" && n_boundary > 0L) {
    cli::cli_abort(c(
      "{n_boundary} indifference point{?s} equal 0 or 1; the linearized transform is undefined there.",
      "i" = "Use {.code boundary = \"clamp\"} (default) or {.code \"drop\"}, or pre-process the data."
    ))
  }
  if (boundary == "clamp") {
    clamped <- d < eps | d > 1 - eps
    d_used <- pmin(pmax(d, eps), 1 - eps)
    n_clamped <- sum(clamped)
  } else if (boundary == "drop") {
    d_used[is_boundary] <- NA_real_
    n_dropped <- n_boundary
  }
  y_lin <- log(1 / d_used - 1) - log(t)
  log_jac <- -log(d_used - d_used^2)          # ln|J|, Eq. 13
  list(
    y_lin = y_lin, d_used = d_used,
    n_boundary = as.integer(n_boundary), n_clamped = as.integer(n_clamped),
    n_dropped = as.integer(n_dropped), log_jac = log_jac
  )
}

#' Per-unit closed-form fit (Hinds et al. 2026, Eqs. 8-9, 12, 14)
#' @keywords internal
#' @noRd
.dd_lin_unit_fit <- function(y_lin, log_jac, conf_level = 0.95) {
  ok <- is.finite(y_lin)
  y <- y_lin[ok]
  lj <- log_jac[ok]
  n <- length(y)
  out <- data.frame(
    n_delays = as.integer(n), logk = NA_real_, se = NA_real_, df = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, k = NA_real_, k_lo = NA_real_, k_hi = NA_real_,
    s2 = NA_real_, loglik_y = NA_real_, loglik_raw = NA_real_
  )
  if (n < 2L) return(out)
  logk <- mean(y)
  rss <- sum((y - logk)^2)
  s2 <- rss / (n - 1)
  se <- sqrt(s2 / n)
  tq <- stats::qt(1 - (1 - conf_level) / 2, df = n - 1)
  sig2_mle <- rss / n
  ll_y <- -(n / 2) * log(2 * pi * sig2_mle) - n / 2
  out$logk <- logk
  out$se <- se
  out$df <- n - 1
  out$ci_lo <- logk - tq * se
  out$ci_hi <- logk + tq * se
  out$k <- exp(logk)
  out$k_lo <- exp(out$ci_lo)
  out$k_hi <- exp(out$ci_hi)
  out$s2 <- s2
  out$loglik_y <- ll_y
  out$loglik_raw <- ll_y + sum(lj)
  out
}

#' One-way random-effects MLEs on the transformed scale (Hinds et al. 2026, Prop. 3.3)
#' @keywords internal
#' @noRd
.dd_lin_re_mle <- function(y_lin, unit, condition, log_jac) {
  unit <- droplevels(as.factor(unit))
  condition <- droplevels(as.factor(condition))
  ok <- is.finite(y_lin)
  y <- y_lin[ok]
  u <- droplevels(unit[ok])
  cnd <- condition[ok]
  lj <- log_jac[ok]
  t_per_unit <- table(u)
  if (length(unique(as.integer(t_per_unit))) != 1L) {
    cli::cli_abort(c(
      "The random-effects model requires a balanced design: every unit must have the same number of delays.",
      "i" = "Counts range {min(t_per_unit)}-{max(t_per_unit)}. Use {.code boundary = \"clamp\"} or complete the data."
    ))
  }
  n_t <- as.integer(t_per_unit[1])
  N <- nlevels(u)
  unit_cond <- factor(tapply(as.character(cnd), u, `[`, 1), levels = levels(cnd))
  ybar_i <- tapply(y, u, mean)
  mu <- tapply(y, cnd, mean)                                       # Eq. 21
  sse_z <- sum((y - ybar_i[as.character(u)])^2)
  ssr_zx <- n_t * sum((ybar_i - mu[as.character(unit_cond)])^2)
  g_zero <- ssr_zx <= 0 || (sse_z / ssr_zx) >= (n_t - 1)
  if (g_zero) {
    sigma2 <- (sse_z + ssr_zx) / (N * n_t)
    g <- 0
  } else {
    sigma2 <- sse_z / (N * (n_t - 1))                              # Eq. 22
    g <- ssr_zx / (sse_z / (n_t - 1)) - 1                          # Eq. 23
  }
  # Eq. 20 at the MLE (third term vanishes at mu-hat), plus the 2*pi constant
  loglik_y <- -(N / 2) * (n_t * log(sigma2) + log1p(g)) -
    (sse_z + ssr_zx / (g + 1)) / (2 * sigma2) - (N * n_t / 2) * log(2 * pi)
  C <- nlevels(cnd)
  list(
    mu = stats::setNames(as.numeric(mu), names(mu)),
    sigma2 = sigma2,
    g = g,
    g_zero = g_zero,
    sigma_u2 = g * sigma2 / n_t,
    sse_z = sse_z,
    ssr_zx = ssr_zx,
    n_units = N,
    n_per_condition = table(unit_cond),
    n_delays = n_t,
    loglik_y = loglik_y,
    loglik_raw = loglik_y + sum(lj),
    n_par = C + 2L,
    unit_mean = stats::setNames(as.numeric(ybar_i), names(ybar_i)),
    unit_condition = unit_cond
  )
}

#' Collapse condition levels declared equal under H0 (Hinds et al. 2026, Sec 3.3 X_red construction)
#' @keywords internal
#' @noRd
.dd_lin_merge_levels <- function(condition, hypothesis = NULL) {
  lv <- levels(condition)
  if (is.null(hypothesis)) hypothesis <- list(lv)
  if (!is.list(hypothesis)) hypothesis <- list(hypothesis)
  bad <- setdiff(unlist(hypothesis), lv)
  if (length(bad)) {
    cli::cli_abort("{.arg hypothesis} names unknown condition level{?s}: {.val {bad}}.")
  }
  if (anyDuplicated(unlist(hypothesis))) {
    cli::cli_abort("A level appears in more than one H0 set.")
  }
  map <- stats::setNames(lv, lv)
  for (set in hypothesis) {
    if (length(set) > 1L) map[set] <- paste(set, collapse = "=")
  }
  factor(map[as.character(condition)], levels = unique(map[lv]))
}

#' Exact F-test for equality of condition means (Hinds et al. 2026, Prop. 3.5)
#'
#' `SSR_{Z|X}` is `T` times the sum over units of the squared difference between
#' the unit mean and its design fitted value. Under `X_full` that fitted value is
#' the condition mean; under `X_red` it is the merged-group mean. The statistic is
#' `F = [(SSE_red - SSE_full) / d_num] / [SSR_{Z|X_full} / d_den]`.
#' @keywords internal
#' @noRd
.dd_lin_ftest <- function(re, hypothesis = NULL) {
  full <- re$unit_condition
  if (nlevels(full) < 2L) {
    cli::cli_abort("No conditions to compare: the fit has a single condition.")
  }
  red <- .dd_lin_merge_levels(full, hypothesis)
  n_t <- re$n_delays
  ybar <- re$unit_mean
  ssr_full <- n_t * sum((ybar - tapply(ybar, full, mean)[as.character(full)])^2)
  ssr_red <- n_t * sum((ybar - tapply(ybar, red, mean)[as.character(red)])^2)
  df1 <- nlevels(full) - nlevels(red)
  df2 <- sum(table(full) - 1)
  f_stat <- ((ssr_red - ssr_full) / df1) / (ssr_full / df2)
  p <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)
  if (re$g_zero) {
    cli::cli_warn("g-hat = 0: Prop. 3.5 assumes g > 0; the F statistic is reported as-is.")
  }
  sets <- if (is.null(hypothesis)) list(levels(full)) else hypothesis
  label <- paste(vapply(sets[lengths(sets) > 1L], paste, "", collapse = " = "), collapse = "; ")
  effect_d <- NA_real_
  if (length(sets) == 1L && length(sets[[1]]) == 2L) {
    m <- re$mu[sets[[1]]]
    effect_d <- sqrt(n_t) * (m[[1]] - m[[2]]) / (sqrt(re$sigma2) * sqrt(re$g + 1))
  }
  tibble::tibble(hypothesis = label, F = f_stat, df1 = df1, df2 = df2, p_value = p, cohens_d = effect_d)
}
