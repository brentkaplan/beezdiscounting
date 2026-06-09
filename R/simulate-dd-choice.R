#' Simulate trial-level SS-vs-LL choice data
#'
#' Generates binary choices from either the structural or descriptive choice
#' model. When `mode = "structural"` (default) choices are drawn from the
#' structural model: `P(LL) = plogis(beta0 + gamma * ((ll/ss) * D(k, delay) -
#' 1))`, `k_i = exp(log_k_pop + sigma_u * u_i)`. When `mode = "descriptive"`
#' choices follow Young's (2018) correlated random-slope logistic model: each
#' subject receives per-subject slopes drawn from a bivariate normal with mean
#' `theta` and covariance `Sigma`, and `P(LL) = plogis((theta[1] + b_i[1]) *
#' log(ll/ss) + (theta[2] + b_i[2]) * log(delay + 1))`. The descriptive branch
#' is the known-truth generator used by recovery tests.
#'
#' @param n_subjects Integer number of subjects.
#' @param ss_amount,ll_amount,delay Numeric trial-design vectors (same length);
#'   each subject is presented all trials. Defaults are a built-in grid.
#' @param mode `"structural"` (default) or `"descriptive"`. Selects the
#'   generative model.
#' @param log_k_pop,sigma_u,gamma,beta0 Structural truth (used only when
#'   `mode = "structural"`).
#' @param equation `"mazur"` or `"exponential"` (structural mode only).
#' @param theta Length-2 numeric vector of population-level slopes for the
#'   descriptive model: `theta[1]` on `log(ll/ss)`, `theta[2]` on
#'   `log(delay + 1)`.
#' @param re_sd Length-2 numeric vector of random-effect standard deviations
#'   (descriptive mode only).
#' @param re_cor Correlation between the two random slopes (descriptive mode
#'   only).
#' @param return_truth Logical. If `TRUE` and `mode = "descriptive"`, the
#'   returned tibble carries two attributes: `subject_slopes` (n_subjects x 2
#'   matrix of realized `b_i` values) and `Sigma` (the 2x2 covariance matrix).
#' @param seed Optional integer seed for reproducibility.
#' @return A [tibble][tibble::tibble] with columns `id`, `ss_amount`,
#'   `ll_amount`, `delay`, `choice` (0/1, 1 = LL chosen). When
#'   `mode = "descriptive"` and `return_truth = TRUE`, the tibble also carries
#'   `subject_slopes` and `Sigma` attributes.
#' @examples
#' # Structural (default)
#' sim <- simulate_dd_choice(n_subjects = 20, seed = 1)
#' head(sim)
#'
#' # Descriptive (Young 2018 correlated random-slope model)
#' sim_desc <- simulate_dd_choice(n_subjects = 20, mode = "descriptive",
#'                                theta = c(1.5, -0.4), seed = 42)
#' head(sim_desc)
#' @export
#' @importFrom stats rnorm rbinom plogis cov
#' @importFrom tibble tibble as_tibble
simulate_dd_choice <- function(n_subjects = 50,
                               ss_amount = c(40, 55, 31, 14, 47, 25, 78, 40, 11, 67),
                               ll_amount = c(65, 75, 85, 25, 60, 30, 80, 65, 30, 75),
                               delay = c(25, 61, 14, 19, 160, 7, 4, 62, 7, 119),
                               mode = c("structural", "descriptive"),
                               log_k_pop = log(0.02), sigma_u = 0.5,
                               gamma = 4, beta0 = 0,
                               equation = c("mazur", "exponential"),
                               theta = c(1.5, -0.4), re_sd = c(0.5, 0.3),
                               re_cor = -0.2, return_truth = FALSE,
                               seed = NULL) {
  mode <- match.arg(mode)
  equation <- match.arg(equation)
  stopifnot(length(ss_amount) == length(ll_amount),
            length(ss_amount) == length(delay))
  if (!is.null(seed)) set.seed(seed)
  n_tr <- length(delay)

  if (mode == "structural") {
    u <- stats::rnorm(n_subjects, 0, sigma_u)
    k <- exp(log_k_pop + u)
    out <- vector("list", n_subjects)
    for (s in seq_len(n_subjects)) {
      D <- if (equation == "mazur") 1 / (1 + k[s] * delay) else exp(-k[s] * delay)
      eta <- beta0 + gamma * ((ll_amount / ss_amount) * D - 1)
      out[[s]] <- tibble::tibble(
        id = as.character(s), ss_amount = ss_amount, ll_amount = ll_amount,
        delay = delay, choice = stats::rbinom(n_tr, 1, stats::plogis(eta)))
    }
    return(tibble::as_tibble(do.call(rbind, out)))
  }

  # descriptive: correlated per-subject slopes on Young's two predictors
  stopifnot(length(theta) == 2L, length(re_sd) == 2L)
  z_mag <- log(ll_amount / ss_amount)
  z_del <- log(delay + 1)
  Sigma <- matrix(c(re_sd[1]^2, re_cor * re_sd[1] * re_sd[2],
                    re_cor * re_sd[1] * re_sd[2], re_sd[2]^2), nrow = 2L)
  Lc <- t(chol(Sigma))
  slopes <- matrix(NA_real_, n_subjects, 2L)
  out <- vector("list", n_subjects)
  for (s in seq_len(n_subjects)) {
    bz <- stats::rnorm(2L)
    bi <- as.numeric(Lc %*% bz)
    slopes[s, ] <- bi
    eta <- (theta[1] + bi[1]) * z_mag + (theta[2] + bi[2]) * z_del
    out[[s]] <- tibble::tibble(
      id = as.character(s), ss_amount = ss_amount, ll_amount = ll_amount,
      delay = delay, choice = stats::rbinom(n_tr, 1, stats::plogis(eta)))
  }
  res <- tibble::as_tibble(do.call(rbind, out))
  if (isTRUE(return_truth)) {
    attr(res, "subject_slopes") <- slopes
    attr(res, "Sigma") <- Sigma
  }
  res
}
