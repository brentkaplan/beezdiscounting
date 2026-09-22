# tests/testthat/test-dd-tmb-guard-info.R
# Guard / clamp / floor activity reporting (audit F-BZ4-2, F-BZ4-5). Reporting
# only: none of these fields changes an estimate.

describe(".dd_discount_mu(clamp = FALSE)", {
  it("returns the raw mean; the default clamps to [1e-6, 1 - 1e-6]", {
    raw <- .dd_discount_mu(0.05, 2920, "exponential", clamp = FALSE)
    expect_equal(raw, exp(-0.05 * 2920))
    expect_lt(raw, 1e-6)
    expect_equal(.dd_discount_mu(0.05, 2920, "exponential"), 1e-6)
  })
})

describe(".dd_tmb_compute_subject_pars latent columns", {
  X <- matrix(1, nrow = 6, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  sid <- rep(0:2, each = 2)
  Sigma <- diag(c(0.25, 4))

  it("phi-target: phi_latent is unfloored, phi is floored at 0.1", {
    co <- c(beta_k = log(0.01), log_phi = log(2))
    b_hat <- cbind(c(0, 0, 0), c(0, -1, -2))   # re_phi = 0, -2, -4
    sp <- .dd_tmb_compute_subject_pars(co, b_hat, c("a", "b", "c"), X, sid,
                                       "mazur", "sltb", n_re = 2L,
                                       Sigma = Sigma, re2_target = 0L)
    expect_named(sp, c("id", "re_k", "re_phi", "k", "phi", "phi_latent"))
    expect_equal(sp$phi_latent, 2 * exp(c(0, -2, -4)))
    expect_equal(sp$phi, pmax(sp$phi_latent, 0.1))
    expect_equal(sum(sp$phi_latent < 0.1), 1L)
  })

  it("s-target: s_latent = exp(log_s + re_s); s is the soft-clamped value", {
    co <- c(beta_k = log(0.01), log_phi = log(10), log_s = 0)
    b_hat <- cbind(c(0, 0, 0), c(0, -1.6, 1.6))  # re_s = 0, -3.2, 3.2
    sp <- .dd_tmb_compute_subject_pars(co, b_hat, c("a", "b", "c"), X, sid,
                                       "green-myerson", "sltb", n_re = 2L,
                                       Sigma = Sigma, re2_target = 1L)
    expect_named(sp, c("id", "re_k", "re_s", "k", "s", "s_latent"))
    expect_equal(sp$s_latent, exp(c(0, -3.2, 3.2)))
    expect_equal(sp$s, .dd_soft_clamp_s_log(log(sp$s_latent)))
    # latent 0.041 / 24.5 sit outside (0.05, 20): clamp active on both.
    act <- .dd_s_clamp_active(sp$s, sp$s_latent)
    expect_identical(act, c(FALSE, TRUE, TRUE))
  })
})

describe(".dd_s_clamp_active", {
  it("flags a latent s AT the lower bound (effective ~0.0518), not interior values", {
    lat <- c(0.05, 1, 20, 0.2)
    eff <- .dd_soft_clamp_s_log(log(lat))
    expect_equal(eff[1], 0.05176, tolerance = 1e-3)
    expect_identical(.dd_s_clamp_active(eff, lat), c(TRUE, FALSE, TRUE, FALSE))
  })
})

describe("fit_dd_tmb guard_info", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("counts rows where the exponential mean guard binds and notes it in summary()", {
    sim <- simulate_dd_ip(n_subjects = 20, equation = "exponential",
                          log_k_pop = log(0.02), sigma_u = 0.5, phi = 12,
                          delays = c(1, 7, 30, 90, 180, 365, 730, 1460),
                          seed = 21)
    fit <- fit_dd_tmb(sim, equation = "exponential", verbose = 0)
    gi <- fit$guard_info
    expect_type(gi, "list")
    expect_gt(gi$mu_guard_lower, 0L)
    expect_identical(gi$n_rows, sum(fit$data$x > 0))
    # independent recount from the fitted subject-level k
    k_row <- .dd_tmb_predict_k(fit, fit$data, level = "subject")
    raw <- exp(-k_row * fit$data$x)
    expect_identical(gi$mu_guard_lower, sum(raw[fit$data$x > 0] < 1e-6))
    notes <- summary(fit)$notes
    expect_true(any(grepl("mean guard", notes)))
  })

  it("reports zero guard activity (and no note) for a typical Mazur fit", {
    sim <- simulate_dd_ip(n_subjects = 20, seed = 1)
    fit <- fit_dd_tmb(sim, equation = "mazur", verbose = 0)
    expect_identical(fit$guard_info$mu_guard_lower, 0L)
    expect_identical(fit$guard_info$mu_guard_upper, 0L)
    expect_null(fit$guard_info$n_s_clamped)
    expect_false(any(grepl("mean guard|clamp|floor", summary(fit)$notes)))
  })

  it("counts clamp-active subjects for k + s ~ 1 and notes them", {
    delays_d <- c(1, 7, 14, 30, 90, 180, 365, 730, 1460, 2920)
    base <- simulate_dd_ip(n_subjects = 18, delays = delays_d,
                           equation = "green-myerson", s = 1.4, sigma_u = 0.5,
                           sigma_s = 0.3, rho_ks = 0.2, phi = 12, seed = 707)
    flat <- data.frame(id = "flat", x = delays_d, y = rep(0.98, length(delays_d)))
    sim <- rbind(data.frame(id = as.character(base$id), x = base$x, y = base$y), flat)
    fit <- fit_dd_tmb(sim, equation = "green-myerson",
                      random_effects = k + s ~ 1, verbose = 0)
    gi <- fit$guard_info
    sp <- fit$subject_pars
    act <- .dd_s_clamp_active(sp$s, sp$s_latent)
    expect_identical(gi$n_s_clamped_lower, sum(act & sp$s < sp$s_latent))
    expect_identical(gi$n_s_clamped_upper, sum(act & sp$s > sp$s_latent))
    # The flat subject sits NEAR the lower bound (latent s ~ 0.06) but the soft
    # clamp moves it by < 1%, so it is not counted.
    expect_lt(sp$s[sp$id == "flat"], 0.1)
    expect_identical(gi$n_s_clamped_lower, 0L)
    expect_false(any(grepl("soft clamp", summary(fit)$notes)))

    # Push one subject's latent s beyond each bound: counted and noted.
    fit2 <- fit
    fit2$subject_pars$s_latent[1:2] <- c(0.03, 30)
    fit2$subject_pars$s[1:2] <- .dd_soft_clamp_s_log(log(c(0.03, 30)))
    gi2 <- .dd_tmb_guard_info(fit2)
    expect_identical(c(gi2$n_s_clamped_lower, gi2$n_s_clamped_upper), c(1L, 1L))
    fit2$guard_info <- gi2
    expect_true(any(grepl("soft clamp on s is active for 2", summary(fit2)$notes)))
  })

  it("stores n_phi_floor for k + phi ~ 1", {
    sim <- simulate_dd_ip(n_subjects = 25, sigma_u = 0.5, sigma_phi = 0.4,
                          rho_kphi = 0.2, phi = 10, seed = 303)
    fit <- fit_dd_tmb(sim, random_effects = k + phi ~ 1, verbose = 0)
    expect_identical(fit$guard_info$n_phi_floor,
                     sum(fit$subject_pars$phi_latent < 0.1))
    fit$subject_pars$phi_latent[1] <- 0.05
    fit$guard_info <- .dd_tmb_guard_info(fit)
    expect_identical(fit$guard_info$n_phi_floor,
                     sum(fit$subject_pars$phi_latent < 0.1))
    expect_gte(fit$guard_info$n_phi_floor, 1L)
    expect_true(any(grepl("phi = 0.1 precision floor", summary(fit)$notes)))
  })
})
