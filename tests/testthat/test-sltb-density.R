# tests/testthat/test-sltb-density.R
# Property tests for the pure-R SLT-beta + Gaussian log-densities.
# These formalize dev/sltb-verification/verify_sltb.R as package tests; the
# C++ template (-nll) is later cross-checked against .dd_slt_logpdf to 1e-8.

S_REF <- 1.0000001
L_REF <- 1e-8

# Analytic truncation normalizer Z (matches verify_sltb.R::Zof)
Zof <- function(mu, phi, s = S_REF, l = L_REF) {
  a <- mu * phi
  b <- (1 - mu) * phi
  pbeta(1 / s + l, a, b) - pbeta(l, a, b)
}

# Exact reference function body from verify_sltb.R (sum of -ll over rows)
ref_nll_mazur <- function(data, par) {
  k <- par[1]
  phi <- par[2]
  delay <- data$delay
  ip <- data$IP
  mu <- 1 / (k * delay + 1)
  alpha <- mu * phi
  beta <- phi * (1 - mu)
  s <- 1.0000001
  l <- 0.00000001
  ll_temp <- lgamma(alpha + beta) - lgamma(beta) - lgamma(alpha) +
    (alpha - 1) * log(ip / s + l) + (beta - 1) * log(1 - (ip / s + l)) -
    log(s) - log(pbeta(1 / s + l, alpha, beta) - pbeta(l, alpha, beta))
  sum(-(ll_temp * ifelse(ip >= 0 & ip <= 1, 1, 0)))
}

describe(".dd_slt_logpdf", {
  mu_grid <- c(0.02, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99)
  phi_grid <- c(1.5, 2, 5, 20, 100, 500)
  grid <- expand.grid(mu = mu_grid, phi = phi_grid)

  it("normalizes: kernel integral equals analytic Z (moderate shapes)", {
    mod <- with(grid, mu * phi >= 1 & (1 - mu) * phi >= 1)
    kern_err <- mapply(function(mu, phi) {
      a <- mu * phi
      b <- (1 - mu) * phi
      num <- stats::integrate(
        function(g) (1 / S_REF) * dbeta(g / S_REF + L_REF, a, b),
        0, 1, rel.tol = 1e-10, subdivisions = 1000L
      )$value
      abs(num - Zof(mu, phi))
    }, grid$mu[mod], grid$phi[mod])
    expect_lt(max(kern_err), 1e-7)
  })

  it("integrates the normalized density to 1 (moderate shapes)", {
    mod <- with(grid, mu * phi >= 1 & (1 - mu) * phi >= 1)
    norm_err <- mapply(function(mu, phi) {
      abs(stats::integrate(
        function(g) exp(.dd_slt_logpdf(g, mu, phi)),
        0, 1, rel.tol = 1e-10
      )$value - 1)
    }, grid$mu[mod], grid$phi[mod])
    expect_lt(max(norm_err), 1e-6)
  })

  it("keeps Z load-bearing (Z varies with shape, not always 1)", {
    expect_lt(Zof(0.05, 1), 0.7)
    expect_gt(Zof(0.5, 10), 0.99)
  })

  it("is finite at the boundaries y=0 and y=1 for all shapes", {
    b0 <- mapply(function(mu, phi) .dd_slt_logpdf(0, mu, phi), grid$mu, grid$phi)
    b1 <- mapply(function(mu, phi) .dd_slt_logpdf(1, mu, phi), grid$mu, grid$phi)
    expect_true(all(is.finite(c(b0, b1))))
  })

  it("matches dbeta/Z in the interior (scale/location limit)", {
    beta_lim_err <- mapply(function(mu, phi) {
      ips <- c(0.1, 0.3, 0.5, 0.7, 0.9)
      max(abs(exp(.dd_slt_logpdf(ips, mu, phi)) -
        dbeta(ips, mu * phi, (1 - mu) * phi) / Zof(mu, phi)))
    }, grid$mu, grid$phi)
    expect_lt(max(beta_lim_err), 1e-3)
  })

  it("equals the reference NLL body exactly (incl. boundary IP=0)", {
    dtest <- data.frame(
      delay = c(7, 30, 180, 365, 730, 1460, 2920),
      IP = c(0.95, 0.8, 0.5, 0.35, 0.2, 0.08, 0.0)
    )
    k0 <- 0.01
    phi0 <- 8
    mu_m <- 1 / (1 + k0 * dtest$delay)
    ours <- sum(-.dd_slt_logpdf(dtest$IP, mu_m, phi0))
    expect_equal(ours, ref_nll_mazur(dtest, c(k0, phi0)), tolerance = 1e-12)
  })

  it("is vectorized over y, mu, and phi", {
    y <- c(0, 0.5, 1)
    mu <- c(0.2, 0.5, 0.8)
    phi <- c(5, 10, 20)
    out <- .dd_slt_logpdf(y, mu, phi)
    expect_length(out, 3L)
    expect_true(all(is.finite(out)))
  })
})

describe(".dd_gaussian_logpdf", {
  it("equals dnorm(log=TRUE)", {
    y <- c(0, 0.25, 0.5, 0.75, 1)
    mu <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    sigma_e <- 0.1
    expect_equal(
      .dd_gaussian_logpdf(y, mu, sigma_e),
      dnorm(y, mu, sigma_e, log = TRUE),
      tolerance = 1e-12
    )
  })

  it("is vectorized and finite for positive sigma_e", {
    out <- .dd_gaussian_logpdf(c(0.2, 0.8), c(0.3, 0.7), 0.05)
    expect_length(out, 2L)
    expect_true(all(is.finite(out)))
  })
})
