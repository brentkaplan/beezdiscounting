# Per-subject two-stage k via NLS (multi-start; single start is fragile).
.two_stage_nls_k <- function(data, equation = "mazur") {
  fo <- if (equation == "mazur") y ~ 1 / (1 + k * x) else y ~ exp(-k * x)
  vapply(split(data, data$id), function(d) {
    for (st in c(exp(-10), 1e-4, 1e-3, 1e-2, 1e-1)) {
      m <- tryCatch(
        stats::nls(fo, data = d, start = list(k = st)),
        error = function(e) NULL
      )
      if (!is.null(m)) return(unname(coef(m)[["k"]]))
    }
    NA_real_
  }, numeric(1))
}

describe("fit_dd_tmb() ties out to per-subject NLS", {
  it("gaussian family agrees with two-stage NLS k on the log scale (mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- .dd_sim_fixture(family = "gaussian", equation = "mazur", seed = 2001)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "gaussian",
                      multi_start = TRUE, verbose = 0)
    sp <- fit$subject_pars                       # per-subject k (BLUP-shrunken)
    ts_raw <- .two_stage_nls_k(sim, equation = "mazur")
    # Align ts to sp$id order: subject_levels uses lex sort ("1","10","11"...)
    # while split() on a factor uses factor level order ("1","2","3"...).
    ts <- ts_raw[sp$id]
    ok <- is.finite(ts) & ts > 0 & is.finite(sp$k) & sp$k > 0
    # gaussian + same mean function => subject k matches NLS k closely in rank/level
    expect_gt(cor(log(sp$k[ok]), log(ts[ok])), 0.95)
  })

  it("population intercept is the GEOMETRIC mean (median), not the arithmetic mean", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- .dd_sim_fixture(family = "sltb", equation = "mazur", seed = 2002)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    k_pop <- exp(unname(fit$model$coefficients[["beta_k"]]))
    ts <- .two_stage_nls_k(sim, equation = "mazur")
    ts <- ts[is.finite(ts) & ts > 0]

    geo  <- exp(mean(log(ts)))     # geometric mean of two-stage k
    med  <- median(ts)            # median (= geometric mean under symmetric log RE)
    arith <- mean(ts)             # arithmetic mean (inflated by exp(sigma^2/2))

    # exp(beta0) tracks the geometric mean / median, within 0.20 relative
    expect_lt(abs(k_pop - geo) / geo, 0.20)
    expect_lt(abs(k_pop - med) / med, 0.25)

    # and it is demonstrably BELOW the arithmetic mean by ~ exp(sigma^2/2):
    # with sigma_u = 0.6, exp(sigma^2/2) = exp(0.18) ~ 1.197, so arith > geo.
    expect_gt(arith, geo)
    sigma_u_hat <- unname(exp(fit$model$coefficients[["log_sigma_u"]]))
    expect_lt(abs((arith / geo) - exp(sigma_u_hat^2 / 2)) / exp(sigma_u_hat^2 / 2), 0.35)
  })
})
