.lin_fit3 <- function() {
  sim <- simulate_dd_linear(n_subjects = 20, delays = c(7, 30, 180, 365, 1095),
                            mu = c(-6.5, -6, -5.5), sigma2 = 2, g = 10, seed = 3)
  fit_dd_linear(sim, factors = "condition")
}
describe("beezdiscounting_linear methods", {
  fit <- .lin_fit3()
  it("print/summary run and mention the clamp bookkeeping", {
    expect_output(print(fit), "Linearized Mazur")
    expect_output(print(summary(fit)), "boundary")
  })
  it("tidy: subject rows by default, population rows with effects='population'", {
    expect_equal(nrow(tidy(fit)), 60L)
    pop <- tidy(fit, effects = "population")
    expect_setequal(pop$term, c("mu_C1", "mu_C2", "mu_C3", "sigma2", "g"))
  })
  it("glance carries mu, sigma2, g, logLik, AIC, BIC with n_par = C + 2", {
    g <- glance(fit)
    expect_equal(g$n_par, 5L)
    expect_equal(g$AIC, -2 * g$logLik + 2 * 5)
    expect_equal(g$BIC, -2 * g$logLik + log(g$nobs) * 5)
  })
  it("coef returns population mu; confint population = ANOVA t-interval with N - C df", {
    expect_equal(coef(fit), fit$re$mu)
    ci <- confint(fit)
    # ORACLE: lm on per-subject ln k computed from the raw data, not from the fit
    d <- fit$data
    logk <- tapply(d$y_lin, d$id, mean)
    cond <- factor(tapply(as.character(d$condition), d$id, `[`, 1), levels = levels(d$condition))
    lm_ci <- confint(lm(logk ~ 0 + cond), level = 0.95)
    expect_equal(unname(ci[, 1]), unname(lm_ci[, 1]), tolerance = 1e-10)
    expect_equal(unname(ci[, 2]), unname(lm_ci[, 2]), tolerance = 1e-10)
  })
  it("augment adds .fitted = 1/(1 + k_i x) per row", {
    a <- augment(fit)
    k <- fit$subjects$k[match(a$id, fit$subjects$id)]
    expect_equal(a$.fitted, 1 / (1 + k * a$x))
  })
  it("logLik defaults: population, raw scale; matches re$loglik_raw; df = C + 2", {
    ll <- logLik(fit)
    expect_equal(as.numeric(ll), fit$re$loglik_raw); expect_equal(attr(ll, "df"), 5L)
    expect_equal(as.numeric(logLik(fit, scale = "transformed")), fit$re$loglik_y)
    expect_equal(as.numeric(logLik(fit, level = "subject")), sum(fit$subjects$loglik_raw))
    expect_equal(attr(logLik(fit, level = "subject"), "df"), 2L * 60L)
  })
  it("nobs is the number of transformed observations", { expect_equal(nobs(fit), 300L) })
  it("anova: default all-equal; pairwise gives 3 rows with cohens_d", {
    a <- anova(fit); expect_equal(nrow(a), 1L); expect_equal(a$df1, 2)
    p <- anova(fit, pairwise = TRUE); expect_equal(nrow(p), 3L); expect_false(anyNA(p$cohens_d))
    expect_s3_class(a, "beezdiscounting_linear_anova")
  })
  it("anova errors when re is NULL or single condition", {
    one <- suppressWarnings(fit_dd_linear(dd_ip))
    expect_error(anova(one), "single condition")
    unb <- suppressWarnings(fit_dd_linear(dd_ip, boundary = "drop"))
    expect_error(anova(unb), "not fitted")
  })
})
