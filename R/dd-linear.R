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
    cli::cli_abort("Delays {.arg x} must be finite and positive; ln(t) is undefined otherwise.")
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
