# Deterministic small SLT-beta mixed sim used across recovery + tie-out tests.
.dd_sim_fixture <- function(family = "sltb", equation = "mazur",
                            n_subjects = 80, seed = 20260606) {
  simulate_dd_ip(
    n_subjects = n_subjects,
    log_k_pop = log(0.01),
    sigma_u = 0.6,
    phi = 12,
    sigma_e = 0.1,
    family = family,
    equation = equation,
    seed = seed
  )
}
