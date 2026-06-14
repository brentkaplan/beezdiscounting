# tests/testthat/test-dd-tmb-phi-re-recovery.R
# 2-RE (log k, log phi) parameter recovery.
#
# Seeds chosen by a documented 12-seed sweep (controller-run; see
# dev/notes/phi-re-sweep.txt). For the design below EVERY swept seed converged
# and recovered all three components inside these tolerances, so the test is not
# seed-fragile; the seeds picked are the most CENTERED in the sampling
# distribution (true sd_u=0.6, sd_phi=0.5, rho=0.4):
#   SEED_SYMM = 10 -> sd_u 0.572, sd_phi 0.484, rho 0.466 (deviations .028/.016/.066)
#   SEED_DIAG =  3 -> sd_u 0.572, sd_phi 0.497           (deviations .028/.003)
# (Policy: never cherry-pick a barely-passing seed.)

describe("2-RE (k, phi) parameter recovery", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  delays <- c(1, 7, 30, 90, 180, 365, 730, 1460, 2920)

  it("recovers sigma_u, sigma_phi, and rho (pdSymm)", {
    sim <- simulate_dd_ip(
      n_subjects = 80, delays = delays,
      sigma_u = 0.6, sigma_phi = 0.5, rho_kphi = 0.4, phi = 12, seed = 10)
    fit <- fit_dd_tmb(sim, random_effects = k + phi ~ 1,
                      covariance_structure = "pdSymm", verbose = 0)
    expect_true(fit$converged)
    vc <- VarCorr(fit)
    expect_equal(vc$StdDev[1], 0.6, tolerance = 0.20)   # sigma_u
    expect_equal(vc$StdDev[2], 0.5, tolerance = 0.25)   # sigma_phi
    expect_equal(vc$Corr[2], 0.4, tolerance = 0.30)     # rho
  })

  it("recovers sigma_u and sigma_phi with rho fixed at 0 (pdDiag)", {
    sim <- simulate_dd_ip(
      n_subjects = 80, delays = delays,
      sigma_u = 0.6, sigma_phi = 0.5, rho_kphi = 0, phi = 12, seed = 3)
    fit <- fit_dd_tmb(sim, random_effects = k + phi ~ 1,
                      covariance_structure = "pdDiag", verbose = 0)
    expect_true(fit$converged)
    vc <- VarCorr(fit)
    expect_equal(vc$StdDev[1], 0.6, tolerance = 0.20)
    expect_equal(vc$StdDev[2], 0.5, tolerance = 0.25)
    expect_true(is.na(vc$Corr[1]) && vc$Corr[2] == 0)   # pdDiag: structural 0
  })

  it("stays finite on a boundary-heavy subject; no phi_i collapses (degeneracy guard)", {
    # A Jarvis-70-style boundary subject (1s then 0s) drives the per-subject SLT
    # likelihood toward phi_i -> 0. The per-subject phi floor (CondExp clamp at
    # 0.1 in the kernel + pmax(., 0.1) in subject_pars) must keep the 2-RE fit
    # finite with every subject's phi_i at or above the floor.
    #
    # NOTE (2026-06-14): unlike the s-target clamp, this one-sided phi floor does NOT
    # exhibit the singular-Hessian HANG, so the smooth-clamp follow-up left it
    # byte-for-byte (still a hard CondExp). Reason: the floor (0.1) is effectively
    # unreachable from data -- three escalating probes (alternating y, then maximally
    # dispersed 0.99/0.01 subjects, then the population precision crushed to ~0.68)
    # bottomed out at min phi_i ~ 0.31, never entering the floor's flat/kinked zone,
    # and every probe converged in 60-81s with no hang. Because the optimum never
    # reaches the kink, the inner-Laplace solve never sees the singular Hessian that
    # hangs the (two-sided, easily-bound) s clamp. See
    # dev/notes/specs/2026-06-14-smooth-clamp-design.md section 4.
    sim <- simulate_dd_ip(n_subjects = 30, sigma_u = 0.5, sigma_phi = 0.4,
                          rho_kphi = 0.2, phi = 10, seed = 303)
    bad <- data.frame(id = "boundary",
                      x = c(7, 30, 180, 365, 730, 1460, 2920),
                      y = c(1,  1,   1,   0,   0,    0,    0))
    sim2 <- rbind(data.frame(id = as.character(sim$id), x = sim$x, y = sim$y), bad)
    fit <- fit_dd_tmb(sim2, random_effects = k + phi ~ 1,
                      covariance_structure = "pdSymm", verbose = 0)
    expect_true(is.finite(fit$loglik))
    expect_true(all(is.finite(fit$subject_pars$phi)))
    expect_true(all(fit$subject_pars$phi >= 0.1 - 1e-6))   # floored, not collapsed
  })
})
