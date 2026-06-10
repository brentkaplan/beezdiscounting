# tests/testthat/test-dd-tmb-1re-regression.R
# Guards the "1-RE path is byte-for-byte unchanged by the 2-RE machinery"
# invariant. The reference coefficients below were captured ONCE on the
# pre-2RE base commit 5f387d2 (develop) for this exact seed/design, via
#   sim <- simulate_dd_ip(n_subjects = 30, sigma_u = 0.6, phi = 10, seed = 99)
#   fit <- fit_dd_tmb(sim, verbose = 0); coef(fit)
# (recorded in dev/notes/ref-1re-5f387d2.txt). simulate_dd_ip()'s default
# sigma_phi = 0 path is RNG-identical to the pre-2RE simulator, so the data is
# the same; if these still match on the 2-RE build, the n_re == 1 kernel/fit
# path is numerically unchanged.

describe("1-RE path is unchanged by the 2-RE machinery", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("default k ~ 1 mazur sltb fit matches its pre-2RE coefficients", {
    sim <- simulate_dd_ip(n_subjects = 30, sigma_u = 0.6, phi = 10, seed = 99)
    fit <- fit_dd_tmb(sim, verbose = 0)
    expect_equal(unname(coef(fit)[["beta_k"]]), -4.6952290113, tolerance = 1e-6)
    expect_equal(unname(coef(fit)[["log_sigma_u"]]), -0.4289307501,
                 tolerance = 1e-6)
    expect_equal(unname(coef(fit)[["log_phi"]]), 2.2816703882, tolerance = 1e-6)
    expect_identical(fit$param_info$n_random_effects, 1L)
  })
})
