describe("fit_dd_choice() structural recovery", {
  it("recovers k, gamma (mazur, intercept off)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_choice(n_subjects = 120, log_k_pop = log(0.02),
                              sigma_u = 0.5, gamma = 4, beta0 = 0,
                              equation = "mazur", seed = 301)
    fit <- fit_dd_choice(sim, mode = "structural", equation = "mazur", verbose = 0)
    expect_true(fit$converged)
    co <- fit$model$coefficients
    k_hat <- exp(unname(co[names(co) == "beta_k"][1]))
    g_hat <- exp(unname(co[["log_gamma"]]))
    expect_lt(abs(k_hat - 0.02) / 0.02, 0.30)
    expect_lt(abs(g_hat - 4) / 4, 0.40)
  })
  it("recovers beta0 when intercept = TRUE", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_choice(n_subjects = 120, log_k_pop = log(0.02),
                              sigma_u = 0.5, gamma = 4, beta0 = 0.6,
                              equation = "mazur", seed = 302)
    fit <- fit_dd_choice(sim, mode = "structural", intercept = TRUE, verbose = 0)
    expect_lt(abs(unname(fit$model$coefficients[["beta0"]]) - 0.6), 0.5)
  })
})

describe("IP <-> choice tie-out", {
  it("structural choice and IP models recover the SAME k (beta0 = 0)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # IP data and choice data from the same k truth
    ip <- simulate_dd_ip(n_subjects = 150, log_k_pop = log(0.02), sigma_u = 0.5,
                         phi = 12, family = "sltb", equation = "mazur", seed = 401)
    ch <- simulate_dd_choice(n_subjects = 150, log_k_pop = log(0.02),
                             sigma_u = 0.5, gamma = 5, beta0 = 0,
                             equation = "mazur", seed = 401)
    fit_ip <- fit_dd_tmb(ip, equation = "mazur", family = "sltb", verbose = 0)
    fit_ch <- fit_dd_choice(ch, mode = "structural", equation = "mazur", verbose = 0)
    k_ip <- exp(unname(fit_ip$model$coefficients[
      names(fit_ip$model$coefficients) == "beta_k"][1]))
    k_ch <- exp(unname(fit_ch$model$coefficients[
      names(fit_ch$model$coefficients) == "beta_k"][1]))
    # both near the truth 0.02, and within a geometric tolerance of each other
    expect_lt(abs(k_ip - 0.02) / 0.02, 0.30)
    expect_lt(abs(k_ch - 0.02) / 0.02, 0.30)
    expect_lt(abs(log(k_ch) - log(k_ip)), log(1.6))   # within ~1.6x of each other
  })
})
