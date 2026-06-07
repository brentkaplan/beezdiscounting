describe(".simulate_dd_ip_mixed()", {
  it("returns the long-format column contract (id, x, y) with no condition by default", {
    sim <- .simulate_dd_ip_mixed(n_subjects = 20, seed = 1)
    expect_s3_class(sim, "tbl_df")
    expect_identical(names(sim), c("id", "x", "y"))
    expect_s3_class(sim$id, "factor")
    expect_equal(nlevels(sim$id), 20L)
    expect_true(all(sim$y >= 0 & sim$y <= 1))
    # 7 default delays per subject
    expect_equal(nrow(sim), 20L * 7L)
  })

  it("adds a condition factor when n_conditions > 1", {
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 12, n_conditions = 2, delta_k = c(0, log(3)), seed = 2
    )
    expect_identical(names(sim), c("id", "condition", "x", "y"))
    expect_s3_class(sim$condition, "factor")
    expect_equal(nlevels(sim$condition), 2L)
  })

  it("errors when delta_k length does not match n_conditions", {
    expect_error(
      .simulate_dd_ip_mixed(n_subjects = 5, n_conditions = 3, delta_k = c(0, 1)),
      "delta_k"
    )
  })

  it("SLT draws recover the population mean curve in expectation (mazur)", {
    set.seed(99)
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 400, log_k_pop = log(0.01), sigma_u = 1e-6, phi = 40,
      family = "sltb", equation = "mazur", seed = 99
    )
    # with sigma_u ~ 0 every subject shares k = exp(log_k_pop); mean y at each delay
    # should track mu = 1/(1 + k*x) to within Monte-Carlo error
    k <- 0.01
    by_delay <- tapply(sim$y, sim$x, mean)
    mu_true <- 1 / (1 + k * as.numeric(names(by_delay)))
    expect_lt(max(abs(by_delay - mu_true)), 0.02)
  })

  it("gaussian draws clamp to [0,1]", {
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 50, family = "gaussian", sigma_e = 0.4, seed = 7
    )
    expect_true(all(sim$y >= 0 & sim$y <= 1))
  })
})

describe(".simulate_dd_ip_mixed() recovery through fit_dd_tmb()", {
  it("recovers population k within 0.15 relative (sltb, mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- .dd_sim_fixture(family = "sltb", equation = "mazur", seed = 101)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    co <- fit$model$coefficients
    k_hat <- exp(unname(co[names(co) == "beta_k"][1]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.15)
  })

  it("recovers population k within 0.30 relative (gaussian, mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # gaussian baseline has [0,1]-clamping bias near the bounds; 0.30 matches P.8
    sim <- .dd_sim_fixture(family = "gaussian", equation = "mazur", seed = 102)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "gaussian",
                      multi_start = TRUE, verbose = 0)
    co <- fit$model$coefficients
    k_hat <- exp(unname(co[names(co) == "beta_k"][1]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.30)
  })

  it("recovers an exponential rate by rank correlation over a k grid", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # exponential is less identifiable; use the identifiable regime + looser bar
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 120, log_k_pop = log(3e-4), sigma_u = 0.6, phi = 12,
      family = "sltb", equation = "exponential", seed = 103
    )
    fit <- fit_dd_tmb(sim, equation = "exponential", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    sp <- fit$subject_pars
    # per-subject two-stage k for the same subjects, ranked agreement
    ts <- vapply(split(sim, sim$id), function(d) {
      m <- tryCatch(stats::nls(y ~ exp(-k * x), data = d, start = list(k = 3e-4)),
                    error = function(e) NULL)
      if (is.null(m)) NA_real_ else unname(coef(m)[["k"]])
    }, numeric(1))
    ts_matched <- ts[as.character(sp$id)]        # align two-stage k to subject_pars order
    ok <- is.finite(ts_matched) & ts_matched > 0
    expect_gt(cor(log(sp$k[ok]), log(ts_matched[ok])), 0.93)
  })
})
