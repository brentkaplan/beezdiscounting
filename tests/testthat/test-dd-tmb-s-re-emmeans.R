# tests/testthat/test-dd-tmb-s-re-emmeans.R
# emmeans beta_k isolation for the (k, s) 2-RE fit.
# NOTE: Do NOT compare a k ~ 1 fit to a k + s ~ 1 fit -- they are different
# likelihoods, so their beta_k need not match (that would be a false invariant).
# Instead assert isolation WITHIN the same s-RE fit: the EMM k_log values equal
# the beta_k linear combinations (so log_sd_re/cor_re did not leak into the
# beta_k vcov slice), and comparisons are finite.
describe("emmeans beta_k isolation under a k + s ~ 1 fit", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  fit2 <- function(seed) {
    sim <- simulate_dd_ip(n_subjects = 40, n_conditions = 2,
                          delta_k = c(0, log(2)),
                          equation = "green-myerson", s = 1.3,
                          sigma_u = 0.5, sigma_s = 0.35, rho_ks = 0.3,
                          phi = 12, seed = seed)
    fit_dd_tmb(sim, equation = "green-myerson",
               random_effects = k + s ~ 1, factors = "condition",
               covariance_structure = "pdSymm", verbose = 0)
  }

  it("get_dd_param_emms isolates beta_k on a (k,s) 2-RE fit", {
    fit <- fit2(7)
    emm <- get_dd_param_emms(fit, factors_in_emm = "condition")
    expect_named(emm, c("level", "k", "k_log", "std.error", "conf.low", "conf.high"))
    # Isolation proof: the EMM k_log values are exactly the beta_k linear
    # combinations (intercept; intercept + conditionC2). If log_sd_re/cor_re had
    # leaked into the beta_k vcov slice, these would not match.
    beta <- unname(fit$model$coefficients[
      names(fit$model$coefficients) == "beta_k"])
    expect_length(beta, 2L)
    expect_equal(sort(emm$k_log), sort(c(beta[1], beta[1] + beta[2])),
                 tolerance = 1e-6)
  })

  it("get_dd_comparisons returns finite contrasts on a (k,s) 2-RE fit", {
    res <- get_dd_comparisons(fit2(8))
    expect_true(all(is.finite(res$k$contrasts_log10$estimate)))
  })
})
