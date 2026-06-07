#' Simulate IP-family mixed-effects discounting data
#'
#' Generates long-format indifference-point data from the mixed-effects
#' discounting model used by [fit_dd_tmb()]. Each subject `i` has a random
#' discount rate `log k_i = log_k_pop + delta_k[condition] + u_i`,
#' `u_i ~ N(0, sigma_u^2)`. The mean indifference proportion at delay `x` is the
#' discounting function `mu` (Mazur hyperbola, exponential, Green-Myerson, or
#' Rachlin), and observed `y`
#' is drawn from the scale-location-truncated beta (`family = "sltb"`) via the
#' inverse-CDF on the truncated beta, or from a clamped Gaussian
#' (`family = "gaussian"`).
#'
#' The SLT draw uses the same constants as the C++ template and the verified
#' reference density: `s_slt = 1.0000001`, `l = 1e-8`, with
#' `a = mu * phi`, `b = (1 - mu) * phi`, and
#' `y = (qbeta(U, a, b) - l) * s_slt` where
#' `U ~ Uniform(pbeta(l, a, b), pbeta(1/s_slt + l, a, b))`.
#'
#' @param n_subjects Integer; number of subjects.
#' @param delays Numeric vector of delays (days) each subject is observed at.
#' @param log_k_pop Numeric; population intercept on `log k`.
#' @param sigma_u Numeric; SD of the subject random intercept on `log k`.
#' @param phi Numeric; SLT-beta precision (`family = "sltb"`).
#' @param sigma_e Numeric; residual SD on `y` (`family = "gaussian"`).
#' @param family One of `"sltb"` (default) or `"gaussian"`.
#' @param equation One of `"mazur"` (default), `"exponential"`,
#'   `"green-myerson"`, or `"rachlin"`. The two 2-parameter forms use the
#'   nonlinearity exponent `s` and reduce to `"mazur"` at `s = 1`.
#' @param s Numeric nonlinearity exponent for the 2-parameter equations
#'   (Green-Myerson / Rachlin); ignored by `"mazur"` / `"exponential"`.
#' @param n_conditions Integer; number of between-subject condition levels. When
#'   `> 1`, a `condition` factor is added and subjects are split across levels.
#' @param delta_k Numeric vector of length `n_conditions`; per-condition shift on
#'   `log k` (the first element is typically `0` for the reference level).
#'   Required (non-`NULL`) when `n_conditions > 1`.
#' @param seed Optional integer seed.
#'
#' @return A [tibble][tibble::tibble] with columns `id` (factor),
#'   `condition` (factor; only when `n_conditions > 1`), `x` (delay), and
#'   `y` (indifference proportion in `[0, 1]`).
#'
#' @examples
#' # A small SLT-beta mixed-effects discounting dataset (id, x, y)
#' sim <- simulate_dd_ip(n_subjects = 8, seed = 1)
#' head(sim)
#'
#' # Two between-subject conditions with a log-k shift on the second level
#' sim2 <- simulate_dd_ip(
#'   n_subjects = 8, n_conditions = 2, delta_k = c(0, log(3)), seed = 2
#' )
#'
#' @export
#' @importFrom stats rnorm runif qbeta pbeta
#' @importFrom tibble tibble
simulate_dd_ip <- function(
  n_subjects = 60,
  delays = c(7, 30, 180, 365, 730, 1460, 2920),
  log_k_pop = log(0.01),
  sigma_u = 0.6,
  phi = 10,
  sigma_e = 0.1,
  family = c("sltb", "gaussian"),
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  s = 1,
  n_conditions = 1,
  delta_k = NULL,
  seed = NULL
) {
  family <- match.arg(family)
  equation <- match.arg(equation)
  if (!is.null(seed)) set.seed(seed)

  if (n_conditions > 1 && is.null(delta_k)) {
    stop("`delta_k` must be supplied (length `n_conditions`) when `n_conditions > 1`.",
         call. = FALSE)
  }
  if (is.null(delta_k)) delta_k <- rep(0, n_conditions)
  if (length(delta_k) != n_conditions) {
    stop("`delta_k` must have length `n_conditions`.", call. = FALSE)
  }

  s_slt <- 1.0000001
  l <- 1e-8

  # assign each subject to a condition (balanced round-robin)
  cond_idx <- rep_len(seq_len(n_conditions), n_subjects)
  cond_lab <- factor(paste0("C", cond_idx), levels = paste0("C", seq_len(n_conditions)))

  u <- stats::rnorm(n_subjects, 0, sigma_u)
  log_k <- log_k_pop + delta_k[cond_idx] + u
  k <- exp(log_k)

  n_delays <- length(delays)
  id <- factor(rep(seq_len(n_subjects), each = n_delays))
  x <- rep(delays, times = n_subjects)
  k_long <- rep(k, each = n_delays)

  mu <- switch(equation,
    mazur = 1 / (1 + k_long * x),
    exponential = exp(-k_long * x),
    `green-myerson` = (1 + k_long * x)^(-s),
    rachlin = ifelse(x > 0, 1 / (1 + k_long * x^s), 1)
  )
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)

  if (family == "sltb") {
    a <- mu * phi
    b <- (1 - mu) * phi
    lo <- stats::pbeta(l, a, b)
    hi <- stats::pbeta(1 / s_slt + l, a, b)
    uu <- stats::runif(length(mu), lo, hi)
    y <- (stats::qbeta(uu, a, b) - l) * s_slt
  } else {
    y <- stats::rnorm(length(mu), mu, sigma_e)
  }
  y <- pmin(pmax(y, 0), 1)

  if (n_conditions > 1) {
    tibble::tibble(
      id = id,
      condition = rep(cond_lab, each = n_delays),
      x = x,
      y = y
    )
  } else {
    tibble::tibble(id = id, x = x, y = y)
  }
}
