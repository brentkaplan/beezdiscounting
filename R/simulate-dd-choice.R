#' Simulate trial-level SS-vs-LL choice data (structural model)
#'
#' Generates binary choices from the structural choice model:
#' `P(LL) = plogis(beta0 + gamma * ((ll/ss) * D(k, delay) - 1))`,
#' `k_i = exp(log_k_pop + sigma_u * u_i)`. Each subject sees the supplied
#' trial design (default = a built-in 10-item SS/LL/delay grid).
#'
#' @param n_subjects Integer.
#' @param ss_amount,ll_amount,delay Numeric trial-design vectors (same length);
#'   each subject is presented all trials. Defaults are a built-in grid.
#' @param log_k_pop,sigma_u,gamma,beta0 Structural truth.
#' @param equation `"mazur"` or `"exponential"`.
#' @param seed Optional integer seed.
#' @return A [tibble][tibble::tibble] with `id`, `ss_amount`, `ll_amount`,
#'   `delay`, `choice` (0/1, 1 = LL).
#' @examples
#' sim <- simulate_dd_choice(n_subjects = 20, seed = 1)
#' head(sim)
#' @export
#' @importFrom stats rnorm rbinom plogis
#' @importFrom tibble tibble as_tibble
simulate_dd_choice <- function(n_subjects = 50,
                               ss_amount = c(40, 55, 31, 14, 47, 25, 78, 40, 11, 67),
                               ll_amount = c(65, 75, 85, 25, 60, 30, 80, 65, 30, 75),
                               delay = c(25, 61, 14, 19, 160, 7, 4, 62, 7, 119),
                               log_k_pop = log(0.02), sigma_u = 0.5,
                               gamma = 4, beta0 = 0,
                               equation = c("mazur", "exponential"),
                               seed = NULL) {
  equation <- match.arg(equation)
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(ss_amount) == length(ll_amount),
            length(ss_amount) == length(delay))
  n_tr <- length(delay)
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
  tibble::as_tibble(do.call(rbind, out))
}
