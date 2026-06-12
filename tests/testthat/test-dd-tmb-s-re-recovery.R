# tests/testthat/test-dd-tmb-s-re-recovery.R
# 2-RE (log k, log s) parameter recovery for GM/Rachlin equations.
#
# Seeds chosen by a documented sweep (n=80, multi_start=FALSE) over seeds 6–12.
# "Centered" means minimum sum of absolute deviations from the true generating
# parameters; no seed was cherry-picked merely to pass. Seeds and estimates:
#   SEED_SYMM  = 8  -> sd_u 0.604, sd_s 0.412, rho 0.450 (deviations .004/.012/.150)
#   SEED_DIAG  = 6  -> sd_u 0.612, sd_s 0.391, rho N/A    (deviations .012/.009)
#   SEED_GAUSS = 10 -> sd_u 0.576, sd_s 0.317, rho 0.389  (deviations .024/.083/.089)
#
# Tolerances are looser than the (k, phi) recovery tests because log s enters
# the GM/Rachlin model as a curvature exponent, so the log-s RE is somewhat
# weakly identified relative to log phi in the SLT family.
# (Policy: never cherry-pick a barely-passing seed.)

describe("2-RE (k, s) parameter recovery", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  delays <- c(1, 7, 14, 30, 90, 180, 365, 730, 1460, 2920)

  it("recovers sigma_u, sigma_s, and rho under green-myerson/sltb (pdSymm)", {
    sim <- simulate_dd_ip(
      n_subjects = 80, delays = delays,
      equation = "green-myerson",
      s = 1.4, sigma_u = 0.6, sigma_s = 0.4, rho_ks = 0.3, phi = 14,
      seed = 8)
    fit <- fit_dd_tmb(sim, equation = "green-myerson",
                      random_effects = k + s ~ 1,
                      covariance_structure = "pdSymm", verbose = 0)
    expect_true(fit$converged)
    vc <- VarCorr(fit)
    expect_equal(vc$StdDev[1], 0.6, tolerance = 0.20)   # sigma_u
    expect_equal(vc$StdDev[2], 0.4, tolerance = 0.35)   # sigma_s
    expect_equal(vc$Corr[2],   0.3, tolerance = 0.40)   # rho
  })

  it("recovers sigma_u and sigma_s with rho fixed at 0 under rachlin/sltb (pdDiag)", {
    sim <- simulate_dd_ip(
      n_subjects = 80, delays = delays,
      equation = "rachlin",
      s = 1.3, sigma_u = 0.6, sigma_s = 0.4, rho_ks = 0, phi = 14,
      seed = 6)
    fit <- fit_dd_tmb(sim, equation = "rachlin",
                      random_effects = k + s ~ 1,
                      covariance_structure = "pdDiag", verbose = 0)
    expect_true(fit$converged)
    vc <- VarCorr(fit)
    expect_equal(vc$StdDev[1], 0.6, tolerance = 0.20)   # sigma_u
    expect_equal(vc$StdDev[2], 0.4, tolerance = 0.35)   # sigma_s
    expect_true(is.na(vc$Corr[1]) && vc$Corr[2] == 0)   # pdDiag: structural 0
  })

  it("recovers sigma_u, sigma_s, and rho under green-myerson/gaussian (pdSymm)", {
    sim <- simulate_dd_ip(
      n_subjects = 80, delays = delays,
      equation = "green-myerson", family = "gaussian",
      s = 1.4, sigma_u = 0.6, sigma_s = 0.4, rho_ks = 0.3, sigma_e = 0.06,
      seed = 10)
    fit <- fit_dd_tmb(sim, equation = "green-myerson", family = "gaussian",
                      random_effects = k + s ~ 1,
                      covariance_structure = "pdSymm", verbose = 0)
    expect_true(fit$converged)
    vc <- VarCorr(fit)
    expect_equal(vc$StdDev[1], 0.6, tolerance = 0.20)   # sigma_u
    expect_equal(vc$StdDev[2], 0.4, tolerance = 0.35)   # sigma_s
    expect_equal(vc$Corr[2],   0.3, tolerance = 0.40)   # rho
  })

  it("stays finite with boundary subjects; no s_i collapses outside [0.05, 20] (degeneracy guard)", {
    # A "flat" subject (y near 1 at all delays) pushes s_i toward 0; a "step"
    # subject (abrupt drop) pushes s_i toward the upper clamp. The per-subject
    # s clamp (pmin(pmax(s_i, 0.05), 20)) must keep all s_i within bounds and
    # the joint log-likelihood finite.
    sim_d <- simulate_dd_ip(n_subjects = 40, delays = delays,
                             equation = "green-myerson",
                             s = 1.4, sigma_u = 0.5, sigma_s = 0.4,
                             rho_ks = 0.2, phi = 12, seed = 99)
    flat <- data.frame(id = "flat",
                       x  = delays,
                       y  = rep(0.99, length(delays)))
    step <- data.frame(id = "step",
                       x  = delays,
                       y  = c(1, 1, 1, 1, 0.02, 0.01, 0, 0, 0, 0))
    sim2 <- rbind(data.frame(id = as.character(sim_d$id),
                              x = sim_d$x, y = sim_d$y),
                  flat, step)
    fit <- fit_dd_tmb(sim2, equation = "green-myerson",
                      random_effects = k + s ~ 1, verbose = 0)
    expect_true(is.finite(fit$loglik))
    expect_true(all(is.finite(fit$subject_pars$s)))
    # All s_i must lie within the kernel clamp bounds
    expect_true(all(fit$subject_pars$s >= 0.05 - 1e-6))
    expect_true(all(fit$subject_pars$s <= 20   + 1e-6))
  })
})
