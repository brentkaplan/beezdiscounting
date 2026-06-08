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
    expect_equal(krow$estimate, exp(unname(fit$model$coefficients["beta_k"][1])),
                 tolerance = 1e-8)
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
})
