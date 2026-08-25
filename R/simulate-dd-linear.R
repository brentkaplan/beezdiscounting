#' Simulate indifference points from the linearized Mazur random-effects model
#'
#' Draws `theta_ic ~ N(mu_c, g sigma2 / T)`, `y_ijc ~ N(theta_ic, sigma2)` and
#' returns `D = 1 / (1 + exp(y) t)` (Hinds et al., 2026, Sec. 4.1), which is
#' always strictly inside (0, 1).
#' @param n_subjects Integer, scalar (recycled) or one per condition.
#' @param delays Positive numeric vector of delays (`T = length(delays)`).
#' @param mu Numeric vector of population mean ln k per condition; names become condition labels.
#' @param sigma2 Measurement-error variance on the transformed scale.
#' @param g Random-effect variance multiplier: `Var(theta) = g * sigma2 / T`.
#' @param seed Optional seed; the global RNG state is restored afterwards.
#' @param attach_truth If `TRUE`, attach `attr(, "truth")`.
#' @return Long `data.frame(id, condition, x, y)`.
#' @examples
#' s <- simulate_dd_linear(30, c(7, 30, 180, 365), mu = c(A = -6, B = -5), sigma2 = 2, g = 10, seed = 1)
#' fit_dd_linear(s, factors = "condition")
#' @export
simulate_dd_linear <- function(n_subjects, delays, mu, sigma2, g, seed = NULL, attach_truth = FALSE) {
  if (any(delays <= 0)) {
    cli::cli_abort("{.arg delays} must be positive.")
  }
  if (sigma2 <= 0 || g < 0) {
    cli::cli_abort("{.arg sigma2} must be > 0 and {.arg g} >= 0.")
  }
  n_cond <- length(mu)
  labels <- if (is.null(names(mu))) paste0("C", seq_len(n_cond)) else names(mu)
  n_subjects <- rep_len(as.integer(n_subjects), n_cond)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
    set.seed(seed)
  }
  n_t <- length(delays)
  cond <- rep(labels, n_subjects)
  ids <- unlist(lapply(seq_len(n_cond), function(c) paste0(labels[c], "_", seq_len(n_subjects[c]))))
  theta <- stats::rnorm(length(ids), mu[match(cond, labels)], sqrt(g * sigma2 / n_t))
  y_lin <- stats::rnorm(length(ids) * n_t, rep(theta, each = n_t), sqrt(sigma2))
  x <- rep(delays, times = length(ids))
  out <- data.frame(
    id = rep(ids, each = n_t), condition = factor(rep(cond, each = n_t), levels = labels),
    x = x, y = 1 / (1 + exp(y_lin) * x), stringsAsFactors = FALSE
  )
  if (attach_truth) {
    attr(out, "truth") <- list(
      theta = stats::setNames(theta, ids), mu = stats::setNames(mu, labels),
      sigma2 = sigma2, g = g, delays = delays
    )
  }
  out
}
