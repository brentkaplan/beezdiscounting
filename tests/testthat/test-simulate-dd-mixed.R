describe("simulate_dd_ip()", {
  it("is an exported function (public simulator API)", {
    expect_true("simulate_dd_ip" %in% getNamespaceExports("beezdiscounting"))
    expect_true(is.function(simulate_dd_ip))
  })

  it("preserves positional arg compatibility for pre-existing formals", {
    # 7th positional must still bind to `family`, not the new sigma_phi
    a <- simulate_dd_ip(8, c(7, 30, 180), log(0.01), 0.6, 10, 0.1, "gaussian",
                        seed = 1)
    b <- simulate_dd_ip(n_subjects = 8, delays = c(7, 30, 180), log_k_pop = log(0.01),
                        sigma_u = 0.6, phi = 10, sigma_e = 0.1, family = "gaussian",
                        seed = 1)
    expect_identical(a, b)
  })

  it("returns the long-format column contract (id, x, y) with no condition by default", {
    sim <- simulate_dd_ip(n_subjects = 20, seed = 1)
    expect_s3_class(sim, "tbl_df")
    expect_identical(names(sim), c("id", "x", "y"))
    expect_s3_class(sim$id, "factor")
    expect_equal(nlevels(sim$id), 20L)
    expect_true(all(sim$y >= 0 & sim$y <= 1))
    # 7 default delays per subject
    expect_equal(nrow(sim), 20L * 7L)
  })

  it("adds a condition factor when n_conditions > 1", {
    sim <- simulate_dd_ip(
      n_subjects = 12, n_conditions = 2, delta_k = c(0, log(3)), seed = 2
    )
    expect_identical(names(sim), c("id", "condition", "x", "y"))
    expect_s3_class(sim$condition, "factor")
    expect_equal(nlevels(sim$condition), 2L)
  })

  it("errors when delta_k length does not match n_conditions", {
    expect_error(
      simulate_dd_ip(n_subjects = 5, n_conditions = 3, delta_k = c(0, 1)),
      "delta_k"
    )
  })

  it("SLT draws recover the population mean curve in expectation (mazur)", {
    set.seed(99)
    sim <- simulate_dd_ip(
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
    sim <- simulate_dd_ip(
      n_subjects = 50, family = "gaussian", sigma_e = 0.4, seed = 7
    )
    expect_true(all(sim$y >= 0 & sim$y <= 1))
  })

  it("green-myerson draws track the (1+k*x)^(-s) mean curve in expectation", {
    set.seed(220)
    s_true <- 0.6
    sim <- simulate_dd_ip(
      n_subjects = 400, log_k_pop = log(0.01), sigma_u = 1e-6, phi = 40,
      family = "sltb", equation = "green-myerson", s = s_true, seed = 220
    )
    k <- 0.01
    by_delay <- tapply(sim$y, sim$x, mean)
    mu_true <- (1 + k * as.numeric(names(by_delay)))^(-s_true)
    expect_lt(max(abs(by_delay - mu_true)), 0.02)
  })

  it("rachlin draws track the 1/(1+k*x^s) mean curve in expectation", {
    set.seed(221)
    s_true <- 1.4
    sim <- simulate_dd_ip(
      n_subjects = 400, log_k_pop = log(0.01), sigma_u = 1e-6, phi = 40,
      family = "sltb", equation = "rachlin", s = s_true, seed = 221
    )
    k <- 0.01
    by_delay <- tapply(sim$y, sim$x, mean)
    mu_true <- 1 / (1 + k * as.numeric(names(by_delay))^s_true)
    expect_lt(max(abs(by_delay - mu_true)), 0.02)
  })

  it("s = 1 makes green-myerson and rachlin identical to mazur (same seed)", {
    a <- simulate_dd_ip(n_subjects = 20, equation = "green-myerson", s = 1,
                        seed = 222)
    b <- simulate_dd_ip(n_subjects = 20, equation = "mazur", seed = 222)
    expect_equal(a$y, b$y, tolerance = 1e-12)
    b2 <- simulate_dd_ip(n_subjects = 20, equation = "rachlin", s = 1, seed = 222)
    expect_equal(b2$y, b$y, tolerance = 1e-12)
  })
})

describe("simulate_dd_ip subject-random phi", {
  it("is unchanged (no phi column) when sigma_phi = 0", {
    a <- simulate_dd_ip(n_subjects = 5, seed = 1)
    b <- simulate_dd_ip(n_subjects = 5, sigma_phi = 0, seed = 1)
    expect_identical(a, b)
    expect_false("phi" %in% names(a))
  })

  it("attaches a per-subject phi column when sigma_phi > 0", {
    sim <- simulate_dd_ip(n_subjects = 8, sigma_phi = 0.5, rho_kphi = 0.3,
                          seed = 1, attach_truth = TRUE)
    expect_true(all(c("id", "x", "y", "phi") %in% names(sim)))
    # phi constant within a subject, varies across subjects
    by_id <- tapply(sim$phi, sim$id, function(v) length(unique(round(v, 10))))
    expect_true(all(by_id == 1L))
    expect_gt(length(unique(round(sim$phi, 6))), 1L)
  })

  it("errors on family = gaussian with sigma_phi > 0", {
    expect_error(simulate_dd_ip(sigma_phi = 0.5, family = "gaussian"),
      regexp = "sltb|gaussian|phi")
  })

  it("does not attach phi when attach_truth = FALSE even if sigma_phi > 0", {
    sim <- simulate_dd_ip(n_subjects = 6, sigma_phi = 0.4, seed = 2)
    expect_false("phi" %in% names(sim))
  })

  it("errors on a non-positive-definite Sigma (sigma_u = 0 or |rho| = 1) when sigma_phi > 0", {
    expect_error(simulate_dd_ip(n_subjects = 6, sigma_u = 0, sigma_phi = 0.4),
                 "positive definite")
    expect_error(simulate_dd_ip(n_subjects = 6, rho_kphi = 1, sigma_phi = 0.4),
                 "positive definite")
  })
})

describe("simulate_dd_ip subject-random s", {
  it("draws a per-subject s when sigma_s > 0 (attach_truth)", {
    sim <- simulate_dd_ip(n_subjects = 40, equation = "green-myerson", s = 1.3,
                          sigma_u = 0.5, sigma_s = 0.3, rho_ks = 0.2,
                          attach_truth = TRUE, seed = 21)
    expect_true("s" %in% names(sim))
    per_subj <- tapply(sim$s, sim$id, function(v) length(unique(round(v, 8))))
    expect_true(all(per_subj == 1))           # constant within subject
    expect_gt(length(unique(round(sim$s, 6))), 1)   # varies across subjects
  })

  it("rejects sigma_s > 0 for a 1-parameter equation", {
    expect_error(
      simulate_dd_ip(n_subjects = 10, equation = "mazur", sigma_s = 0.3, seed = 1),
      "green-myerson|rachlin")
  })

  it("rejects sigma_s and sigma_phi both > 0 (q = 2 only)", {
    expect_error(
      simulate_dd_ip(n_subjects = 10, equation = "green-myerson",
                     sigma_s = 0.3, sigma_phi = 0.3, seed = 1),
      "both")
  })

  it("appends sigma_s/rho_ks after attach_truth (positional API stability)", {
    # New args MUST be appended at the end so every pre-existing positional slot
    # (incl. attach_truth) keeps its index. Assert the formals order directly --
    # robust, unlike a brittle 15-positional call.
    fa <- names(formals(simulate_dd_ip))
    expect_lt(which(fa == "attach_truth"), which(fa == "sigma_s"))
    expect_lt(which(fa == "attach_truth"), which(fa == "rho_ks"))
    expect_equal(tail(fa, 2), c("sigma_s", "rho_ks"))
  })

  it("does not attach s when attach_truth = FALSE even if sigma_s > 0", {
    sim <- simulate_dd_ip(n_subjects = 6, equation = "green-myerson",
                          sigma_s = 0.4, seed = 3)
    expect_false("s" %in% names(sim))
  })

  it("errors on non-positive-definite Sigma when sigma_s > 0", {
    expect_error(simulate_dd_ip(n_subjects = 6, equation = "green-myerson",
                                sigma_u = 0, sigma_s = 0.4), "positive definite")
    expect_error(simulate_dd_ip(n_subjects = 6, equation = "green-myerson",
                                rho_ks = 1, sigma_s = 0.4), "positive definite")
  })

  it("draws a per-subject s under rachlin too", {
    sim <- simulate_dd_ip(n_subjects = 30, equation = "rachlin", s = 1.3,
                          sigma_u = 0.5, sigma_s = 0.3, attach_truth = TRUE, seed = 5)
    expect_true("s" %in% names(sim))
    per_subj <- tapply(sim$s, sim$id, function(v) length(unique(round(v, 8))))
    expect_true(all(per_subj == 1))
    expect_gt(length(unique(round(sim$s, 6))), 1)
  })
})

describe("simulate_dd_ip() recovery through fit_dd_tmb()", {
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
    sim <- simulate_dd_ip(
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
