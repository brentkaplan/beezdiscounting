.skip_unless_full_tests()
# tests/testthat/test-dd-tmb-phi-re-recovery.R
# 2-RE (log k, log phi) parameter recovery.
#
# Protocol (audit F-BZ4-6, 2026-09-22): same as test-dd-tmb-s-re-recovery.R.
# Prespecified seeds 101:108 (frozen before running; seeds 1-12 had been swept
# earlier and are not used); every seed must converge with a PD Hessian;
# recovery is judged on the mean over seeds against tolerances fixed before the
# run; failures are reported by seed, never re-selected. The boundary-subject
# regression lives in test-dd-tmb-2re-boundary.R (ungated).

recovery_seeds <- 101:108

run_phi_recovery <- function(sim_args, cov) {
  res <- lapply(recovery_seeds, function(sd) {
    sim <- do.call(simulate_dd_ip, c(sim_args, list(seed = sd)))
    fit <- fit_dd_tmb(sim, random_effects = k + phi ~ 1,
                      covariance_structure = cov, verbose = 0)
    vc <- VarCorr(fit)
    data.frame(
      seed = sd,
      converged = isTRUE(fit$converged),
      hessian_pd = isTRUE(fit$hessian_pd),
      sd_u = vc$StdDev[1], sd_phi = vc$StdDev[2], rho = vc$Corr[2]
    )
  })
  do.call(rbind, res)
}

describe("2-RE (k, phi) parameter recovery (prespecified seeds 101:108)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  delays <- c(1, 7, 30, 90, 180, 365, 730, 1460, 2920)

  it("recovers sigma_u, sigma_phi, and rho (pdSymm)", {
    tab <- run_phi_recovery(
      list(n_subjects = 80, delays = delays, sigma_u = 0.6, sigma_phi = 0.5,
           rho_kphi = 0.4, phi = 12), "pdSymm")
    bad <- tab$seed[!(tab$converged & tab$hessian_pd)]
    expect_true(length(bad) == 0L,
                info = paste("not converged / non-PD Hessian at seeds:",
                             paste(bad, collapse = ", ")))
    info <- paste(utils::capture.output(print(tab)), collapse = "\n")
    expect_equal(mean(tab$sd_u),   0.6, tolerance = 0.20, info = info)
    expect_equal(mean(tab$sd_phi), 0.5, tolerance = 0.25, info = info)
    expect_equal(mean(tab$rho),    0.4, tolerance = 0.30, info = info)
  })

  it("recovers sigma_u and sigma_phi with rho fixed at 0 (pdDiag)", {
    tab <- run_phi_recovery(
      list(n_subjects = 80, delays = delays, sigma_u = 0.6, sigma_phi = 0.5,
           rho_kphi = 0, phi = 12), "pdDiag")
    bad <- tab$seed[!(tab$converged & tab$hessian_pd)]
    expect_true(length(bad) == 0L,
                info = paste("not converged / non-PD Hessian at seeds:",
                             paste(bad, collapse = ", ")))
    info <- paste(utils::capture.output(print(tab)), collapse = "\n")
    expect_equal(mean(tab$sd_u),   0.6, tolerance = 0.20, info = info)
    expect_equal(mean(tab$sd_phi), 0.5, tolerance = 0.25, info = info)
    expect_true(all(tab$rho == 0))   # pdDiag: structural 0
  })
})
