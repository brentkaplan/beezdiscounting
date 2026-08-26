describe(".dd_lin_unit_fit()", {
  y <- c(-2.1, -2.6, -2.35, -2.9)          # already-transformed y_j (Eq. 6)
  lj <- c(1.2, 1.5, 1.3, 2.0)              # ln|J_j|
  it("ln k is the sample mean of y (Eq. 8) and k = exp(ln k) (Eq. 9)", {
    r <- .dd_lin_unit_fit(y, lj)
    expect_equal(r$logk, mean(y))
    expect_equal(r$k, exp(mean(y)))
    expect_equal(r$n_delays, 4L)
  })
  it("S^2 uses divisor T-1, SE = S/sqrt(T), CI uses t_{T-1} (Sec 2.2)", {
    r <- .dd_lin_unit_fit(y, lj, conf_level = 0.95)
    s2 <- sum((y - mean(y))^2) / 3
    expect_equal(r$s2, s2)
    expect_equal(r$se, sqrt(s2 / 4))
    expect_equal(r$df, 3)
    tq <- qt(0.975, df = 3)
    expect_equal(r$ci_lo, mean(y) - tq * sqrt(s2 / 4))
    expect_equal(r$ci_hi, mean(y) + tq * sqrt(s2 / 4))
    expect_equal(r$k_lo, exp(r$ci_lo)); expect_equal(r$k_hi, exp(r$ci_hi))
  })
  it("log-liks: transformed scale uses MLE variance (divisor T); raw adds sum ln|J| (Eq. 12-14)", {
    r <- .dd_lin_unit_fit(y, lj)
    sig2_mle <- sum((y - mean(y))^2) / 4
    ll_y <- sum(dnorm(y, mean(y), sqrt(sig2_mle), log = TRUE))
    expect_equal(r$loglik_y, ll_y)
    expect_equal(r$loglik_raw, ll_y + sum(lj))
  })
  it("drops NA y (drop-mode boundary points) before estimating", {
    r <- .dd_lin_unit_fit(c(y, NA), c(lj, NA))
    expect_equal(r$logk, mean(y)); expect_equal(r$n_delays, 4L)
  })
  it("returns NA row with n_delays when fewer than 2 usable points", {
    r <- .dd_lin_unit_fit(c(-2, NA), c(1, NA))
    expect_equal(r$n_delays, 1L); expect_true(is.na(r$logk)); expect_true(is.na(r$se))
  })
  it("geometric-mean identity: k = geomean((1/D - 1)/t)", {
    d <- c(0.9, 0.6, 0.3); t <- c(1, 10, 100)
    tr <- .dd_lin_transform(d, t)
    r <- .dd_lin_unit_fit(tr$y_lin, tr$log_jac)
    expect_equal(r$k, exp(mean(log((1 / d - 1) / t))))
  })
})
