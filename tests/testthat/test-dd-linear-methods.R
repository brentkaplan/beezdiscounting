.lin_fit3 <- function() {
  sim <- simulate_dd_linear(n_subjects = 20, delays = c(7, 30, 180, 365, 1095),
                            mu = c(-6.5, -6, -5.5), sigma2 = 2, g = 10, seed = 3)
  fit_dd_linear(sim, factors = "condition")
}
.lin_fit_ncz0 <- function() {
  d <- data.frame(
    id = rep(c("a", "b"), each = 3), x = rep(c(1, 7, 30), 2),
    y = c(0.9, 0.5, 0.2, 0.8, 0.4, 0.1), cond = rep(c("A", "B"), each = 3)
  )
  fit_dd_linear(d, factors = "cond")
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
    # ORACLE: lm on per-subject ln k computed from the fit's stored y_lin (independent of the
    # RE closed forms, not of the transform)
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
  it("confint population errors when N - C = 0 (one subject per condition)", {
    fit0 <- .lin_fit_ncz0()
    expect_error(confint(fit0), "N - C")
  })
  it("anova errors when the F-test denominator df is 0 (N - C = 0)", {
    fit0 <- .lin_fit_ncz0()
    expect_error(anova(fit0), "denominator df")
  })
})

describe("predict() and .std_resid for beezdiscounting_linear", {
  fit <- .lin_fit3()
  it("predict(type = 'parameters') renames the subject table to the plot contract", {
    sp <- predict(fit, type = "parameters")
    expect_s3_class(sp, "tbl_df")
    expect_equal(
      names(sp),
      c("id", "condition", "logk", "logk_lower", "logk_upper", "k", "k_lower", "k_upper")
    )
    expect_equal(sp$k, fit$subjects$k)
    expect_equal(sp$k_lower, fit$subjects$k_lo)
    expect_equal(sp$k_upper, fit$subjects$k_hi)
    expect_equal(sp$logk_lower, fit$subjects$ci_lo)
    expect_identical(predict(fit), sp)
  })
  it("predict rejects newdata and other types", {
    expect_error(predict(fit, newdata = fit$data[1:3, ]), "augment")
    expect_error(predict(fit, type = "response"), "should be")
  })
  it("augment .std_resid = (y_lin - logk_i) / sqrt(re$sigma2)", {
    a <- augment(fit)
    logk <- fit$subjects$logk[match(a$id, fit$subjects$id)]
    expect_equal(a$.std_resid, (a$y_lin - logk) / sqrt(fit$re$sigma2))
    expect_equal(as.vector(tapply(a$.std_resid, a$id, mean)), rep(0, 60), tolerance = 1e-10)
  })
  it("augment .std_resid uses the pooled within-subject MS when re is NULL; NA at dropped rows", {
    unb <- suppressWarnings(fit_dd_linear(dd_ip, boundary = "drop"))
    a <- augment(unb)
    expect_true(all(is.na(a$.std_resid[!is.finite(a$y_lin)])))
    s <- unb$subjects
    pooled <- sum(s$df * s$s2) / sum(s$df)
    logk <- s$logk[match(a$id, s$id)]
    ok <- is.finite(a$y_lin)
    expect_equal(a$.std_resid[ok], ((a$y_lin - logk) / sqrt(pooled))[ok])
  })
  it("augment .std_resid is not identically +/- 1/sqrt(2) for two-delay subjects", {
    two <- fit_dd_linear(simulate_dd_linear(
      n_subjects = 6, delays = c(7, 30), mu = -5, sigma2 = 2, g = 5, seed = 4
    ))
    r <- abs(augment(two)$.std_resid)
    expect_false(all(abs(r - 1 / sqrt(2)) < 1e-8))
  })
})

describe("confint() default level: 0.95 for population, conf_level for subject", {
  sim <- simulate_dd_linear(n_subjects = 10, delays = c(7, 30, 180, 365), mu = c(-6.5, -6),
                            sigma2 = 2, g = 10, seed = 5)
  fit90 <- fit_dd_linear(sim, factors = "condition", conf_level = 0.90)
  it("parm = 'subject' reproduces subjects$ci_lo / ci_hi", {
    ci <- confint(fit90, parm = "subject")
    expect_equal(unname(ci[, 1]), fit90$subjects$ci_lo)
    expect_equal(unname(ci[, 2]), fit90$subjects$ci_hi)
  })
  it("population interval defaults to 0.95 like the other tiers; explicit level still wins", {
    expect_equal(confint(fit90), confint(fit90, level = 0.95))
    expect_false(isTRUE(all.equal(confint(fit90), confint(fit90, level = 0.90))))
    expect_equal(confint(fit90, parm = "subject", level = 0.95)[1, ],
                 confint(fit90, level = 0.95, parm = "subject")[1, ])
    expect_false(isTRUE(all.equal(
      unname(confint(fit90, parm = "subject", level = 0.95)[, 1]), fit90$subjects$ci_lo
    )))
  })
})
