.skip_unless_full_tests()
# tests/testthat/test-dd-tmb-s-re-recovery.R
# 2-RE (log k, log s) parameter recovery for GM/Rachlin equations.
#
# Protocol (audit F-BZ4-6, 2026-09-22). Each condition is fit on the SAME
# prespecified seeds, 101:108, frozen before any of them was run (seeds 1-17 had
# been examined by earlier sweeps and are not used). Per seed the fit must be
# converged with a positive-definite Hessian and at most one subject may have an
# active s soft clamp (.dd_s_clamp_active(): effective vs latent s differ by
# > 1% on the log scale). Recovery is judged on the MEAN over seeds against the
# tolerances below, which were fixed before the run. A failing seed or condition
# is reported by seed number, never re-selected; a failure is a finding to
# investigate, not a reason to change seeds or tolerances.
#
# Tolerances are looser than the (k, phi) recovery tests because log s enters the
# GM/Rachlin model as a curvature exponent, so the log-s RE is weakly identified.
# Heavy (8 fits per condition): runs only under BEEZ_FULL_TESTS=true.

recovery_seeds <- 101:108

run_s_recovery <- function(sim_args, fit_args) {
  res <- lapply(recovery_seeds, function(sd) {
    sim <- do.call(simulate_dd_ip, c(sim_args, list(seed = sd)))
    fit <- do.call(fit_dd_tmb, c(list(sim), fit_args, list(verbose = 0)))
    sp <- fit$subject_pars
    vc <- VarCorr(fit)
    data.frame(
      seed = sd,
      converged = isTRUE(fit$converged),
      hessian_pd = isTRUE(fit$hessian_pd),
      n_clamped = sum(.dd_s_clamp_active(sp$s, sp$s_latent)),
      sd_u = vc$StdDev[1], sd_s = vc$StdDev[2], rho = vc$Corr[2]
    )
  })
  do.call(rbind, res)
}

expect_recovery_protocol <- function(tab) {
  bad_conv <- tab$seed[!(tab$converged & tab$hessian_pd)]
  expect_true(length(bad_conv) == 0L,
              info = paste("not converged / non-PD Hessian at seeds:",
                           paste(bad_conv, collapse = ", ")))
  bad_clamp <- tab$seed[tab$n_clamped > 1L]
  expect_true(length(bad_clamp) == 0L,
              info = paste("more than one clamp-active subject at seeds:",
                           paste(bad_clamp, collapse = ", ")))
}

describe("2-RE (k, s) parameter recovery (prespecified seeds 101:108)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  delays <- c(1, 7, 14, 30, 90, 180, 365, 730, 1460, 2920)

  it("recovers sigma_u, sigma_s, and rho under green-myerson/sltb (pdSymm)", {
    tab <- run_s_recovery(
      list(n_subjects = 80, delays = delays, equation = "green-myerson",
           s = 1.4, sigma_u = 0.6, sigma_s = 0.4, rho_ks = 0.3, phi = 14),
      list(equation = "green-myerson", random_effects = k + s ~ 1,
           covariance_structure = "pdSymm"))
    expect_recovery_protocol(tab)
    info <- paste(utils::capture.output(print(tab)), collapse = "\n")
    expect_equal(mean(tab$sd_u), 0.6, tolerance = 0.20, info = info)
    expect_equal(mean(tab$sd_s), 0.4, tolerance = 0.35, info = info)
    expect_equal(mean(tab$rho),  0.3, tolerance = 0.40, info = info)
  })

  it("recovers sigma_u and sigma_s with rho fixed at 0 under rachlin/sltb (pdDiag)", {
    tab <- run_s_recovery(
      list(n_subjects = 80, delays = delays, equation = "rachlin",
           s = 1.3, sigma_u = 0.6, sigma_s = 0.4, rho_ks = 0, phi = 14),
      list(equation = "rachlin", random_effects = k + s ~ 1,
           covariance_structure = "pdDiag"))
    expect_recovery_protocol(tab)
    info <- paste(utils::capture.output(print(tab)), collapse = "\n")
    expect_equal(mean(tab$sd_u), 0.6, tolerance = 0.20, info = info)
    expect_equal(mean(tab$sd_s), 0.4, tolerance = 0.35, info = info)
    expect_true(all(tab$rho == 0))   # pdDiag: structural 0
  })

  it("recovers sigma_u, sigma_s, and rho under green-myerson/gaussian (pdSymm)", {
    tab <- run_s_recovery(
      list(n_subjects = 80, delays = delays, equation = "green-myerson",
           family = "gaussian", s = 1.4, sigma_u = 0.6, sigma_s = 0.4,
           rho_ks = 0.3, sigma_e = 0.06),
      list(equation = "green-myerson", family = "gaussian",
           random_effects = k + s ~ 1, covariance_structure = "pdSymm"))
    expect_recovery_protocol(tab)
    info <- paste(utils::capture.output(print(tab)), collapse = "\n")
    expect_equal(mean(tab$sd_u), 0.6, tolerance = 0.20, info = info)
    expect_equal(mean(tab$sd_s), 0.4, tolerance = 0.35, info = info)
    expect_equal(mean(tab$rho),  0.3, tolerance = 0.40, info = info)
  })
})
