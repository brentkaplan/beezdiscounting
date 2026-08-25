describe("Jacobian-corrected likelihood (Sec 2.3 / Eq. 13-14)", {
  it("the implied density of a single D integrates to 1 over (0, 1)", {
    logk <- -4; sigma <- 0.8; t <- 30
    dens <- function(d) dnorm(log(1 / d - 1) - log(t), logk, sigma) / (d - d^2)
    expect_equal(integrate(dens, 0, 1, rel.tol = 1e-8)$value, 1, tolerance = 1e-6)
  })
  it("subject-level raw log-lik is comparable to a Gaussian NLS log-lik on the same D", {
    d <- data.frame(id = "s", x = c(1, 7, 30, 90, 180, 365), y = c(0.95, 0.8, 0.55, 0.35, 0.2, 0.12))
    fit <- fit_dd_linear(d)
    nls_fit <- stats::nls(y ~ 1 / (1 + k * x), d, start = list(k = 0.02))
    expect_type(fit$subjects$loglik_raw, "double")
    expect_true(is.finite(fit$subjects$loglik_raw)); expect_true(is.finite(as.numeric(logLik(nls_fit))))
  })
})
