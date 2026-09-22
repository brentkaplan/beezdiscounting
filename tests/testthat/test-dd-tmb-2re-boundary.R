# tests/testthat/test-dd-tmb-2re-boundary.R
# Boundary-subject regressions for the 2-RE fits (moved out of the gated
# recovery files so they run in every devtools::test()).

describe("2-RE boundary subjects", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("(k, s): fits finitely, converged with a PD Hessian, on a near-bound subject", {
    # Regression for the singular-Hessian HANG. A "no discounting" subject (flat y
    # across every delay, including 2920) drives its s_i toward the LOWER bound
    # (s -> 0). Under the old HARD clamp that landed in the zero-gradient/kinked
    # zone -> singular inner-Laplace Hessian -> the fit ground to iter_max. The
    # SOFT clamp (.dd_soft_clamp_s_log / kernel logspace_add) removes the kink.
    # The "step" subject adds a sharply discounting subject to stress the joint
    # fit. Audit F-BZ4-6: the fit must also be converged with a PD Hessian.
    delays_d <- c(1, 7, 14, 30, 90, 180, 365, 730, 1460, 2920)
    base <- simulate_dd_ip(n_subjects = 18, delays = delays_d,
                           equation = "green-myerson",
                           s = 1.4, sigma_u = 0.5, sigma_s = 0.3,
                           rho_ks = 0.2, phi = 12, seed = 707)
    flat <- data.frame(id = "flat", x = delays_d, y = rep(0.98, length(delays_d)))
    step <- data.frame(id = "step", x = delays_d,
                       y = c(0.99, 0.99, 0.99, 0.99, 0.5, 0.02, 0.02, 0.02, 0.02, 0.02))
    sim <- rbind(data.frame(id = as.character(base$id), x = base$x, y = base$y),
                 flat, step)
    fit <- fit_dd_tmb(sim, equation = "green-myerson",
                      random_effects = k + s ~ 1, verbose = 0)
    expect_true(fit$converged)
    expect_true(fit$hessian_pd)
    expect_true(is.finite(fit$loglik))
    expect_true(all(is.finite(fit$subject_pars$s)))
    expect_true(all(fit$subject_pars$s >= 0.05 - 1e-9 &
                    fit$subject_pars$s <= 20   + 1e-9))
    # The flat subject is driven NEAR the lower bound (s < 0.1). Its latent
    # value (~0.06) stays inside the bound, so the soft clamp barely moves it
    # (guard_info counts it only when the move exceeds 1% on the log scale).
    expect_lt(fit$subject_pars$s[fit$subject_pars$id == "flat"], 0.1)
  })

  it("(k, phi): stays finite on a boundary-heavy subject; no phi_i collapses", {
    # A Jarvis-70-style boundary subject (1s then 0s) drives the per-subject SLT
    # likelihood toward phi_i -> 0. The per-subject phi floor (CondExp clamp at
    # 0.1 in the kernel + pmax(., 0.1) in subject_pars) must keep the 2-RE fit
    # finite with every subject's phi_i at or above the floor.
    #
    # NOTE (2026-06-14): unlike the s-target clamp, this one-sided phi floor does NOT
    # exhibit the singular-Hessian HANG, so the smooth-clamp follow-up left it
    # byte-for-byte (still a hard CondExp). The floor (0.1) is effectively
    # unreachable from data -- escalating probes bottomed out at min phi_i ~ 0.31.
    # See dev/notes/specs/2026-06-14-smooth-clamp-design.md section 4.
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
    expect_identical(fit$guard_info$n_phi_floor, 0L)
  })
})
