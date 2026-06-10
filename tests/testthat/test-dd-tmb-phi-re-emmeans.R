# tests/testthat/test-dd-tmb-phi-re-emmeans.R
describe("emmeans beta_k isolation under a 2-RE fit", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  fit2 <- function(seed) {
    sim <- simulate_dd_ip(n_subjects = 40, n_conditions = 2,
                          delta_k = c(0, log(2)), sigma_u = 0.5,
                          sigma_phi = 0.4, rho_kphi = 0.3, phi = 12, seed = seed)
    fit_dd_tmb(sim, random_effects = k + phi ~ 1, factors = "condition",
               covariance_structure = "pdSymm", verbose = 0)
  }

  it("get_dd_param_emms isolates beta_k (tibble) on a 2-RE fit", {
    fit <- fit2(7)
    emm <- get_dd_param_emms(fit, factors_in_emm = "condition")
    expect_s3_class(emm, "tbl_df")
    expect_named(emm, c("level", "k", "k_log", "std.error",
                        "conf.low", "conf.high"))
    expect_equal(nrow(emm), 2L)
    expect_true(all(is.finite(emm$k)))
    # Isolation proof: the EMM k_log values are exactly the beta_k linear
    # combinations (intercept; intercept + conditionC2). If log_sd_re/cor_re
    # had leaked into the beta slice, these would not match.
    beta <- unname(fit$model$coefficients[
      names(fit$model$coefficients) == "beta_k"])
    expect_length(beta, 2L)
    expect_equal(sort(emm$k_log), sort(c(beta[1], beta[1] + beta[2])),
                 tolerance = 1e-6)
  })

  it("get_dd_comparisons returns finite contrasts on a 2-RE fit", {
    res <- get_dd_comparisons(fit2(8))   # auto-detects the condition factor
    expect_s3_class(res, "beezdiscounting_comparison")
    expect_named(res$k, c("emmeans", "contrasts_log10", "contrasts_ratio"))
    expect_true(all(is.finite(res$k$contrasts_log10$estimate)))
  })
})
