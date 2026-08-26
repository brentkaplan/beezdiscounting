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
    # ORACLE: computed only from the raw `d` frame, never from fit$data or .dd_lin_transform.
    yv <- log(1 / d$y - 1) - log(d$x)
    s2 <- mean((yv - mean(yv))^2)
    ll_raw_oracle <- sum(dnorm(yv, mean(yv), sqrt(s2), log = TRUE)) + sum(-log(d$y - d$y^2))
    expect_equal(fit$subjects$loglik_raw, ll_raw_oracle, tolerance = 1e-10)
    expect_equal(fit$subjects$loglik_raw - fit$subjects$loglik_y, sum(-log(d$y - d$y^2)), tolerance = 1e-10)
  })
})
