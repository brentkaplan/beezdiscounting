.lin_re_sim <- function(n = c(30, 40), T = 5, mu = c(-6, -5), sigma2 = 2, g = 8, seed = 1) {
  set.seed(seed)
  cond <- rep(c("A", "B"), n)
  theta <- stats::rnorm(sum(n), mu[match(cond, c("A", "B"))], sqrt(g * sigma2 / T))
  data.frame(
    unit = factor(rep(seq_len(sum(n)), each = T)),
    condition = factor(rep(cond, each = T)),
    y_lin = stats::rnorm(sum(n) * T, rep(theta, each = T), sqrt(sigma2)),
    log_jac = stats::runif(sum(n) * T, 1, 3)
  )
}

describe(".dd_lin_re_mle()", {
  d <- .lin_re_sim()
  fit <- .dd_lin_re_mle(d$y_lin, d$unit, d$condition, d$log_jac)

  it("mu_c is the condition mean of all transformed obs (Eq. 21)", {
    expect_equal(fit$mu, tapply(d$y_lin, d$condition, mean)[c("A", "B")], ignore_attr = TRUE)
    expect_named(fit$mu, c("A", "B"))
  })
  it("SSE_Z and SSR_Z|X are the within-unit and T x between-unit sums of squares", {
    ybar_i <- tapply(d$y_lin, d$unit, mean)
    ybar_c <- tapply(d$y_lin, d$condition, mean)
    cond_i <- tapply(as.character(d$condition), d$unit, `[`, 1)
    expect_equal(fit$sse_z, sum((d$y_lin - ybar_i[as.character(d$unit)])^2))
    expect_equal(fit$ssr_zx, 5 * sum((ybar_i - ybar_c[cond_i])^2))
  })
  it("sigma2 and g follow Eqs. 22-23 when SSE_Z/SSR < T-1", {
    N <- 70; T <- 5
    expect_false(fit$g_zero)
    expect_equal(fit$sigma2, fit$sse_z / (N * (T - 1)))
    expect_equal(fit$g, fit$ssr_zx / (fit$sse_z / (T - 1)) - 1)
    expect_equal(fit$sigma_u2, fit$g * fit$sigma2 / T)
    expect_equal(fit$n_par, 4L)
  })
  it("ORACLE: matches nlme::lme(method = 'ML') on the same transformed data", {
    skip_on_cran()
    lme_fit <- nlme::lme(y_lin ~ 0 + condition, random = ~ 1 | unit, data = d, method = "ML")
    expect_equal(unname(nlme::fixef(lme_fit)), unname(fit$mu), tolerance = 1e-6)
    expect_equal(lme_fit$sigma^2, fit$sigma2, tolerance = 1e-4)
    vc <- as.numeric(nlme::VarCorr(lme_fit)[1, "Variance"])
    expect_equal(vc, fit$sigma_u2, tolerance = 1e-3)
    expect_equal(as.numeric(logLik(lme_fit)), fit$loglik_y, tolerance = 1e-6)
  })
  it("g = 0 branch: sigma2 = (SSE_Z + SSR)/(N T) when SSE_Z/SSR >= T-1", {
    set.seed(2)
    dd <- .lin_re_sim(n = c(20, 20), T = 4, mu = c(0, 0), sigma2 = 1, g = 0)
    # force no between-unit spread beyond noise: shrink unit means toward condition means
    ybar_i <- tapply(dd$y_lin, dd$unit, mean); ybar_c <- tapply(dd$y_lin, dd$condition, mean)
    ci <- tapply(as.character(dd$condition), dd$unit, `[`, 1)
    dd$y_lin <- dd$y_lin - ybar_i[as.character(dd$unit)] + ybar_c[ci[as.character(dd$unit)]] +
      0.01 * (ybar_i[as.character(dd$unit)] - ybar_c[ci[as.character(dd$unit)]])
    f0 <- .dd_lin_re_mle(dd$y_lin, dd$unit, dd$condition, dd$log_jac)
    expect_true(f0$g_zero); expect_equal(f0$g, 0)
    # ORACLE (independent): at g = 0 the model is iid N(mu_c, sigma2) over all N*T obs,
    # i.e. lm(y ~ 0 + condition) with ML variance (logLik.lm uses the ML divisor).
    lm0 <- stats::lm(y_lin ~ 0 + condition, data = dd)
    expect_equal(f0$sigma2, sum(stats::resid(lm0)^2) / (40 * 4))
    expect_equal(f0$loglik_y, as.numeric(stats::logLik(lm0)), tolerance = 1e-8)
    expect_equal(unname(f0$mu), unname(stats::coef(lm0)), tolerance = 1e-10)
  })
  it("raw-scale log-lik = transformed log-lik + sum ln|J| (Sec 3.4)", {
    expect_equal(fit$loglik_raw, fit$loglik_y + sum(d$log_jac))
  })
  it("aborts when every unit's transformed points are identical (SSE_Z = 0)", {
    dz <- data.frame(
      unit = factor(rep(seq_len(4), each = 4)),
      condition = factor(rep(c("A", "B"), each = 8)),
      y_lin = rep(c(-6, -5.5, -5, -6.2), each = 4),
      log_jac = 1
    )
    expect_error(.dd_lin_re_mle(dz$y_lin, dz$unit, dz$condition, dz$log_jac), "SSE_Z = 0")
  })
  it("aborts on unbalanced T", {
    d2 <- d[-1, ]
    expect_error(.dd_lin_re_mle(d2$y_lin, d2$unit, d2$condition, d2$log_jac), "balanced")
  })
})
