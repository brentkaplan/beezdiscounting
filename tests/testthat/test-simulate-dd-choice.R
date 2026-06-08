describe("simulate_dd_choice() structural", {
  it("is exported and returns the generic choice frame", {
    expect_true("simulate_dd_choice" %in% getNamespaceExports("beezdiscounting"))
    sim <- simulate_dd_choice(n_subjects = 10, seed = 1)
    expect_s3_class(sim, "tbl_df")
    expect_true(all(c("id", "ss_amount", "ll_amount", "delay", "choice") %in% names(sim)))
    expect_true(all(sim$choice %in% c(0, 1)))
  })
  it("P(LL) tracks the structural curve in expectation (sigma_u ~ 0)", {
    set.seed(7)
    sim <- simulate_dd_choice(n_subjects = 400, log_k_pop = log(0.02),
                              sigma_u = 1e-6, gamma = 4, beta0 = 0,
                              equation = "mazur", seed = 7)
    k <- 0.02; gamma <- 4
    D <- 1 / (1 + k * sim$delay)
    p_true <- plogis(gamma * ((sim$ll_amount / sim$ss_amount) * D - 1))
    by_bin <- tapply(sim$choice, cut(p_true, 5), mean)
    centers <- tapply(p_true, cut(p_true, 5), mean)
    expect_lt(max(abs(by_bin - centers), na.rm = TRUE), 0.07)
  })
})
