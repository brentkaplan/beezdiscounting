skip_on_cran(); skip_if_not_installed("TMB")
.ccm <- new.env(parent = emptyenv())
.get_choice_fit <- function() {
  if (!exists("fit", envir = .ccm)) {
    # reuse the Task-4 fixture builder (helper-dd-choice.R)
    dat <- .choice_fit_fixture(n_subjects = 40, gamma = 4, seed = 21)
    .ccm$fit <- fit_dd_choice(dat, mode = "structural", equation = "mazur",
                              verbose = 0)
  }
  .ccm$fit
}

describe("beezdiscounting_choice S3", {
  it("coef/fixef return the raw named optimizer vector incl. log_gamma", {
    fit <- .get_choice_fit()
    co <- coef(fit)
    expect_true(all(c("beta_k", "log_sigma_u", "log_gamma") %in% names(co)))
    expect_identical(nlme::fixef(fit), co)
  })
  it("ranef returns id/u_i/k", {
    fit <- .get_choice_fit()
    re <- ranef(fit)
    expect_named(re, c("id", "u_i", "k"))
    expect_equal(nrow(re), fit$param_info$n_subjects)
  })
  it("predict(type='prob') returns P(LL) in (0,1) of length nobs", {
    fit <- .get_choice_fit()
    pr <- predict(fit, type = "prob")
    expect_true(".prob" %in% names(pr) || "predict.id" %in% names(pr))
    p <- if (".prob" %in% names(pr)) pr$.prob else pr$predict.id
    expect_length(p, nrow(fit$data)); expect_true(all(p > 0 & p < 1))
  })
  it("predict(type='parameters') returns the subject k table", {
    fit <- .get_choice_fit()
    pp <- predict(fit, type = "parameters")
    expect_true(all(c("id", "u_i", "k") %in% names(pp)))
  })
  it("tidy returns the 8-col broom contract with a k:(Intercept) row", {
    fit <- .get_choice_fit()
    td <- tidy(fit)
    expect_named(td, c("term", "estimate", "std.error", "statistic",
                       "p.value", "component", "estimate_scale", "term_display"))
    expect_true("k:(Intercept)" %in% td$term)
  })
  it("tidy natural space exponentiates k and reports gamma on natural scale", {
    fit <- .get_choice_fit()
    tn <- tidy(fit, report_space = "natural")
    krow <- tn[tn$term == "k:(Intercept)", ]
    beta_k1 <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1])
    expect_equal(krow$estimate, exp(beta_k1), tolerance = 1e-8)
    grow <- tn[tn$term == "gamma", ]
    expect_equal(grow$estimate, exp(unname(fit$model$coefficients[["log_gamma"]])),
                 tolerance = 1e-8)
  })
  it("glance reports mode/equation/backend and one row", {
    fit <- .get_choice_fit()
    g <- glance(fit)
    expect_equal(nrow(g), 1L)
    expect_equal(g$mode, "structural"); expect_equal(g$equation, "mazur")
  })
  it("confint natural exponentiates the k row; summary/print run", {
    fit <- .get_choice_fit()
    ci <- confint(fit, parm = "k:(Intercept)", report_space = "natural")
    expect_true(ci$conf.low <= ci$estimate && ci$estimate <= ci$conf.high)
    expect_s3_class(summary(fit), "summary.beezdiscounting_choice")
    expect_invisible(print(fit))
  })
  it("logLik/AIC/BIC/nobs match stored values", {
    fit <- .get_choice_fit()
    expect_equal(as.numeric(logLik(fit)), fit$loglik)
    expect_equal(AIC(fit), fit$AIC); expect_equal(BIC(fit), fit$BIC)
    expect_equal(nobs(fit), fit$param_info$n_obs)
  })
  it("augment appends finite .fitted/.resid/.std_resid", {
    fit <- .get_choice_fit()
    au <- augment(fit)
    expect_true(all(c(".fitted", ".resid", ".std_resid") %in% names(au)))
    expect_equal(nrow(au), nrow(fit$data))
    expect_true(all(is.finite(au$.fitted)))
    expect_true(all(is.finite(au$.resid)))
    expect_true(all(is.finite(au$.std_resid)))
  })
  it("tidy reports a beta0 shape row on the identity scale when intercept=TRUE", {
    if (!exists("fit_on", envir = .ccm)) {
      dat <- .choice_fit_fixture(n_subjects = 40, gamma = 4, beta0 = 0.5, seed = 22)
      .ccm$fit_on <- fit_dd_choice(dat, mode = "structural", equation = "mazur",
                                   intercept = TRUE, verbose = 0)
    }
    fit_on <- .ccm$fit_on
    td_log <- tidy(fit_on)                       # internal/default space
    td_nat <- tidy(fit_on, report_space = "natural")
    brow_log <- td_log[td_log$term == "beta0", ]
    brow_nat <- td_nat[td_nat$term == "beta0", ]
    expect_equal(nrow(brow_log), 1L)
    expect_equal(brow_log$component, "shape")
    expect_equal(brow_log$estimate_scale, "identity")
    expect_equal(brow_log$estimate, unname(fit_on$model$coefficients[["beta0"]]),
                 tolerance = 1e-8)
    # identity scale: estimate must NOT change across report spaces
    expect_equal(brow_nat$estimate, brow_log$estimate, tolerance = 1e-12)
  })
})
