# Shared deterministic structural-choice fixture for the choice family tests.
.choice_fit_fixture <- function(n_subjects = 40, log_k_pop = log(0.02),
                                sigma_u = 0.5, gamma = 4, beta0 = 0,
                                equation = "mazur", seed = 1) {
  set.seed(seed)
  ss_grid <- c(40, 55, 31, 14, 47, 25, 78, 40, 11, 67)
  ll_grid <- c(65, 75, 85, 25, 60, 30, 80, 65, 30, 75)
  delay_grid <- c(25, 61, 14, 19, 160, 7, 4, 62, 7, 119)
  u <- rnorm(n_subjects, 0, sigma_u)
  k <- exp(log_k_pop + u)
  rows <- do.call(rbind, lapply(seq_len(n_subjects), function(s) {
    D <- if (equation == "mazur") 1 / (1 + k[s] * delay_grid) else exp(-k[s] * delay_grid)
    eta <- beta0 + gamma * ((ll_grid / ss_grid) * D - 1)
    data.frame(id = as.character(s), ss_amount = ss_grid, ll_amount = ll_grid,
               delay = delay_grid, choice = rbinom(length(eta), 1, plogis(eta)))
  }))
  rows
}

# Deterministic descriptive-choice fixture (correlated random slopes).
.choice_desc_fixture <- function(n_subjects = 80, theta = c(1.6, -0.45),
                                 re_sd = c(0.5, 0.3), re_cor = -0.3, seed = 11) {
  simulate_dd_choice(
    n_subjects = n_subjects, mode = "descriptive",
    theta = theta, re_sd = re_sd, re_cor = re_cor, seed = seed)
}
