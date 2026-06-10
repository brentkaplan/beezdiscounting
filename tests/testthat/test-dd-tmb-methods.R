# tests/testthat/test-dd-tmb-methods.R
skip_on_cran()
skip_if_not_installed("TMB")

# ---------------------------------------------------------------------------
# File-level memoized fit: compile + fit ONCE, reuse across all tests.
# Mirror beezdemand's .fdt_cache pattern.
# ---------------------------------------------------------------------------
.ddm_cache <- new.env(parent = emptyenv())

.get_fit_for_methods <- function() {
  if (!exists("fit", envir = .ddm_cache)) {
    dat <- .dd_sim_fixture(family = "sltb", equation = "mazur", n_subjects = 30,
                           seed = 101)
    .ddm_cache$fit <- fit_dd_tmb(
      dat,
      equation       = "mazur",
      family         = "sltb",
      random_effects = k ~ 1,
      multi_start    = TRUE,
      verbose        = 0
    )
  }
  .ddm_cache$fit
}

.get_gm_fit_for_methods <- function() {
  if (!exists("gm_fit", envir = .ddm_cache)) {
    dat <- simulate_dd_ip(n_subjects = 40,
                          delays = c(1, 7, 30, 180, 365, 730, 1460, 2920),
                          log_k_pop = log(0.01), sigma_u = 0.5, phi = 15,
                          s = 0.6, family = "sltb",
                          equation = "green-myerson", seed = 501)
    .ddm_cache$gm_fit <- fit_dd_tmb(dat, equation = "green-myerson",
                                    family = "sltb", multi_start = TRUE,
                                    verbose = 0)
  }
  .ddm_cache$gm_fit
}

# ---------------------------------------------------------------------------
# M.0: shared helpers
# ---------------------------------------------------------------------------

describe(".dd_tmb_build_term_names()", {
  it("maps beta_k columns to k:<design colname> and passes aux through", {
    fit <- .get_fit_for_methods()
    tn <- .dd_tmb_build_term_names(fit)
    expect_true("k:(Intercept)" %in% tn$term)
    # aux + sigma_u rows keep their raw display names
    expect_true(any(tn$term %in% c("log_phi", "log_sigma_e")))
    expect_true("log_sigma_u" %in% tn$term)
    expect_length(tn$k_idx, ncol(fit$formula_details$X))
    expect_type(tn$other_idx, "integer")
  })

  it("maps log_s to the display term s for a 2-parameter fit", {
    fit <- .get_gm_fit_for_methods()
    tn <- .dd_tmb_build_term_names(fit)
    expect_true("s" %in% tn$term)
    expect_false("log_s" %in% tn$term)   # raw name is mapped to display "s"
  })
})

describe(".dd_tmb_model_se()", {
  it("returns a named SE vector aligned to coefficients", {
    fit <- .get_fit_for_methods()
    se <- .dd_tmb_model_se(fit)
    co <- fit$model$coefficients
    expect_named(se, names(co))
    expect_length(se, length(co))
  })
})

# ---------------------------------------------------------------------------
# M.1: logLik / AIC / BIC / nobs
# ---------------------------------------------------------------------------

describe("logLik / AIC / BIC / nobs", {
  it("logLik carries df and nobs attributes", {
    fit <- .get_fit_for_methods()
    ll <- logLik(fit)
    expect_s3_class(ll, "logLik")
    expect_equal(as.numeric(ll), fit$loglik)
    expect_equal(attr(ll, "df"), length(fit$opt$par))
    expect_equal(attr(ll, "nobs"), fit$param_info$n_obs)
  })

  it("AIC matches the stored value and honours custom k", {
    fit <- .get_fit_for_methods()
    expect_equal(AIC(fit), fit$AIC)
    np <- length(fit$opt$par)
    expect_equal(AIC(fit, k = 3), 2 * (-fit$loglik) + 3 * np)
  })

  it("BIC and nobs match stored values", {
    fit <- .get_fit_for_methods()
    expect_equal(BIC(fit), fit$BIC)
    expect_equal(nobs(fit), fit$param_info$n_obs)
  })
})

# ---------------------------------------------------------------------------
# M.2: coef / fixef / ranef
# ---------------------------------------------------------------------------

describe("coef / fixef", {
  it("coef() returns the raw named optimizer vector", {
    fit <- .get_fit_for_methods()
    co <- coef(fit)
    expect_type(co, "double")
    expect_named(co)
    expect_true("beta_k" %in% names(co))
    expect_true("log_sigma_u" %in% names(co))
    expect_true(any(c("log_phi", "log_sigma_e") %in% names(co)))
    expect_identical(co, fit$model$coefficients)
  })

  it("fixef() returns the same named numeric vector as coef()", {
    fit <- .get_fit_for_methods()
    expect_identical(nlme::fixef(fit), coef(fit))
  })

  it("coef() returns log_s for a 2-parameter (green-myerson) fit", {
    fit <- .get_gm_fit_for_methods()
    expect_true("log_s" %in% names(coef(fit)))
  })
})

describe("ranef", {
  it("returns id, u_i, and per-subject k for every subject", {
    fit <- .get_fit_for_methods()
    re <- ranef(fit)
    expect_s3_class(re, "data.frame")
    expect_named(re, c("id", "u_i", "k"))
    expect_equal(nrow(re), fit$param_info$n_subjects)
    expect_equal(re$id, fit$subject_pars$id)
    expect_equal(re$u_i, fit$subject_pars$u_i)
    expect_true(all(re$k > 0))
  })

  it("omits phi: phi is population-level, not a subject-level column", {
    fit <- .get_fit_for_methods()
    re <- ranef(fit)
    expect_false("phi" %in% names(re))
    expect_named(re, c("id", "u_i", "k"))
  })
})

# ---------------------------------------------------------------------------
# M.3: predict
# ---------------------------------------------------------------------------

describe("predict", {
  it("type = 'parameters' returns the subject_pars tibble", {
    fit <- .get_fit_for_methods()
    pp <- predict(fit, type = "parameters")
    expect_s3_class(pp, "tbl_df")
    expect_true(all(c("id", "u_i", "k") %in% names(pp)))
    expect_equal(nrow(pp), fit$param_info$n_subjects)
  })

  it("type = 'response', level = 'subject' adds .fitted in (0,1)", {
    fit <- .get_fit_for_methods()
    pr <- predict(fit, type = "response", level = "subject")
    expect_true(".fitted" %in% names(pr))
    expect_equal(nrow(pr), nrow(fit$data))
    expect_true(all(pr$.fitted > 0 & pr$.fitted < 1))
  })

  it("level = 'subject' equals the discounting fn at the subject's k", {
    fit <- .get_fit_for_methods()
    pr <- predict(fit, type = "response", level = "subject")
    sp <- fit$subject_pars
    k_by_id <- stats::setNames(sp$k, sp$id)
    x <- fit$data[[fit$param_info$x_var]]
    id <- as.character(fit$data[[fit$param_info$id_var]])
    mu_raw <- 1 / (1 + k_by_id[id] * x)  # mazur
    mu_guarded <- pmin(pmax(mu_raw, 1e-6), 1 - 1e-6)
    expect_equal(pr$.fitted, unname(mu_guarded), tolerance = 1e-8)
  })

  it("level = 'population' uses RE = 0 and needs no id column", {
    fit <- .get_fit_for_methods()
    nd <- data.frame(x = c(7, 30, 180, 365))
    pr <- predict(fit, newdata = nd, type = "response", level = "population")
    expect_true("predict.fixed" %in% names(pr))
    k_pop <- exp(unname(fit$model$coefficients["beta_k"][1]))
    mu_raw <- 1 / (1 + k_pop * nd$x)
    mu_guarded <- pmin(pmax(mu_raw, 1e-6), 1 - 1e-6)
    expect_equal(pr$predict.fixed, mu_guarded, tolerance = 1e-8)
  })

  it("population predictions differ from subject predictions when RE nonzero", {
    fit <- .get_fit_for_methods()
    sigma_u <- exp(fit$model$coefficients[["log_sigma_u"]])
    # RE variance is nonzero → subject and population predictions must differ
    # on the full dataset
    if (sigma_u > 1e-4) {
      pr_s <- predict(fit, type = "response", level = "subject")
      pr_p <- predict(fit, type = "response", level = "population")
      expect_false(isTRUE(all.equal(pr_s$.fitted, pr_p$predict.fixed)))
    }
  })

  it("level = c('population','subject') returns both columns", {
    fit <- .get_fit_for_methods()
    nd <- fit$data
    pr <- predict(fit, newdata = nd, type = "response",
                  level = c("population", "subject"))
    expect_true(all(c("predict.fixed", "predict.id") %in% names(pr)))
  })

  it("rejects a numeric nlme-style level", {
    fit <- .get_fit_for_methods()
    expect_error(predict(fit, level = 1), "should be one of")
  })

  it("exponential fit yields .fitted = exp(-k * x) at population level", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    d2 <- simulate_dd_ip(n_subjects = 25, family = "gaussian",
                                equation = "exponential", seed = 7)
    f2 <- fit_dd_tmb(d2, equation = "exponential", family = "gaussian",
                     random_effects = k ~ 1, verbose = 0)
    pr <- predict(f2, type = "response", level = "population",
                  newdata = data.frame(x = c(7, 365)))
    k_pop <- exp(unname(f2$model$coefficients["beta_k"][1]))
    mu_raw <- exp(-k_pop * c(7, 365))
    mu_guarded <- pmin(pmax(mu_raw, 1e-6), 1 - 1e-6)
    expect_equal(pr$predict.fixed, mu_guarded, tolerance = 1e-8)
  })

  it("ERRORS on an unseen factor level in newdata (no silent zero-pad)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    d3 <- simulate_dd_ip(n_subjects = 30, family = "sltb",
                                equation = "mazur", seed = 31)
    d3$grp <- factor(rep(c("ctrl", "trt"),
                         length.out = length(unique(d3$id)))[
                           match(d3$id, unique(d3$id))])
    f3 <- fit_dd_tmb(d3, factors = "grp", verbose = 0)
    nd <- data.frame(x = c(7, 30),
                     id = f3$param_info$subject_levels[1],
                     grp = factor("NEVER_SEEN"))
    expect_error(
      predict(f3, newdata = nd, type = "response", level = "population"),
      regexp = "not seen in the fit|unseen"
    )
  })

  it("green-myerson predictions use (1+k*x)^(-s) at population level", {
    fit <- .get_gm_fit_for_methods()
    nd <- data.frame(x = c(7, 180, 730))
    pr <- predict(fit, newdata = nd, type = "response", level = "population")
    k_pop <- exp(unname(fit$model$coefficients[
      names(fit$model$coefficients) == "beta_k"][1]))
    s_hat <- exp(unname(fit$model$coefficients[["log_s"]]))
    mu_raw <- (1 + k_pop * nd$x)^(-s_hat)
    expect_equal(pr$predict.fixed, pmin(pmax(mu_raw, 1e-6), 1 - 1e-6),
                 tolerance = 1e-8)
  })
})

# ---------------------------------------------------------------------------
# M.4: augment / fitted / residuals
# ---------------------------------------------------------------------------

describe("augment", {
  it("adds .fitted, .resid, .std_resid and preserves row count", {
    fit <- .get_fit_for_methods()
    aug <- augment(fit)
    expect_s3_class(aug, "tbl_df")
    expect_true(all(c(".fitted", ".resid", ".std_resid") %in% names(aug)))
    expect_equal(nrow(aug), nrow(fit$data))
    y <- fit$data[[fit$param_info$y_var]]
    expect_equal(aug$.resid, y - aug$.fitted, tolerance = 1e-10)
  })

  it(".std_resid is finite and positive-signed (resid / sd > 0 iff resid > 0)", {
    fit <- .get_fit_for_methods()
    aug <- augment(fit)
    expect_true(all(is.finite(aug$.std_resid)))
    # sign of .std_resid matches sign of .resid
    expect_true(all(sign(aug$.resid) == sign(aug$.std_resid) |
                      aug$.resid == 0))
  })

  it("supports newdata and computes resid against the y_var in newdata", {
    fit <- .get_fit_for_methods()
    nd <- fit$data[seq_len(10L), , drop = FALSE]
    aug <- augment(fit, newdata = nd)
    expect_equal(nrow(aug), 10L)
    expect_true(all(c(".fitted", ".resid") %in% names(aug)))
  })

  it("gaussian fit uses constant sigma_e for .std_resid", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    d2 <- simulate_dd_ip(n_subjects = 25, family = "gaussian",
                                equation = "mazur", seed = 9)
    f2 <- fit_dd_tmb(d2, equation = "mazur", family = "gaussian",
                     random_effects = k ~ 1, verbose = 0)
    aug <- augment(f2)
    sigma_e <- exp(f2$model$coefficients[["log_sigma_e"]])
    expect_equal(aug$.std_resid, aug$.resid / sigma_e, tolerance = 1e-10)
  })
})

describe("fitted / residuals", {
  it("fitted() returns a numeric vector of length nobs", {
    fit <- .get_fit_for_methods()
    fv <- stats::fitted(fit)
    expect_type(fv, "double")
    expect_length(fv, nobs(fit))
    expect_true(all(fv > 0 & fv < 1))
  })

  it("residuals(type = 'response') equals y - fitted", {
    fit <- .get_fit_for_methods()
    y <- fit$data[[fit$param_info$y_var]]
    r <- stats::residuals(fit, type = "response")
    expect_equal(r, y - stats::fitted(fit), tolerance = 1e-10)
  })

  it("residuals(type = 'pearson') are finite", {
    fit <- .get_fit_for_methods()
    rp <- stats::residuals(fit, type = "pearson")
    expect_true(all(is.finite(rp)))
  })
})

# ---------------------------------------------------------------------------
# M.5: tidy
# ---------------------------------------------------------------------------

describe("tidy", {
  it("returns exactly the 8-column broom contract in spec order", {
    fit <- .get_fit_for_methods()
    td <- tidy(fit)
    expect_s3_class(td, "tbl_df")
    expect_named(td, c("term", "estimate", "std.error", "statistic",
                       "p.value", "component", "estimate_scale",
                       "term_display"))
  })

  it("has nrow > 0 and no NA estimates on a good fit (fixed rows)", {
    fit <- .get_fit_for_methods()
    td <- tidy(fit, effects = "fixed")
    expect_gt(nrow(td), 0L)
    expect_true(all(!is.na(td$estimate)))
  })

  it("fixed-effect rows carry component = 'fixed', variance rows 'variance'", {
    fit <- .get_fit_for_methods()
    td <- tidy(fit)
    expect_true(any(td$component == "fixed"))
    expect_true(any(td$component == "variance"))
    k_row <- td[td$term == "k:(Intercept)", ]
    expect_equal(nrow(k_row), 1L)
    expect_equal(k_row$component, "fixed")
  })

  it("internal space leaves the log-k estimate untransformed", {
    fit <- .get_fit_for_methods()
    td <- tidy(fit, report_space = "internal")
    k_row <- td[td$term == "k:(Intercept)", ]
    expect_equal(k_row$estimate, unname(fit$model$coefficients["beta_k"][1]),
                 tolerance = 1e-10)
    expect_equal(k_row$estimate_scale, "log")
  })

  it("natural space exponentiates the intercept to k and rescales SE", {
    fit <- .get_fit_for_methods()
    td_int <- tidy(fit, report_space = "internal")
    td_nat <- tidy(fit, report_space = "natural")
    ki <- td_int[td_int$term == "k:(Intercept)", ]
    kn <- td_nat[td_nat$term == "k:(Intercept)", ]
    expect_equal(kn$estimate, exp(ki$estimate), tolerance = 1e-8)
    expect_equal(kn$std.error, exp(ki$estimate) * ki$std.error, tolerance = 1e-8)
    expect_equal(kn$estimate_scale, "natural")
  })

  it("non-circular: natural-scale k equals exp(beta_k) computed independently", {
    fit <- .get_fit_for_methods()
    td_nat <- tidy(fit, report_space = "natural")
    kn <- td_nat[td_nat$term == "k:(Intercept)", ]
    expected_k <- exp(unname(fit$model$coefficients["beta_k"][1]))
    expect_equal(kn$estimate, expected_k, tolerance = 1e-8)
  })

  it("log10 space divides the log estimate by log(10)", {
    fit <- .get_fit_for_methods()
    td_int <- tidy(fit, report_space = "internal")
    td_l10 <- tidy(fit, report_space = "log10")
    ki <- td_int[td_int$term == "k:(Intercept)", ]
    kl <- td_l10[td_l10$term == "k:(Intercept)", ]
    expect_equal(kl$estimate, ki$estimate / log(10), tolerance = 1e-8)
    expect_equal(kl$estimate_scale, "log10")
  })

  it("statistic/p.value are estimation-scale invariant to report_space", {
    fit <- .get_fit_for_methods()
    a <- tidy(fit, report_space = "internal")
    b <- tidy(fit, report_space = "natural")
    c_ <- tidy(fit, report_space = "log10")
    ka <- a[a$term == "k:(Intercept)", ]
    kb <- b[b$term == "k:(Intercept)", ]
    kc <- c_[c_$term == "k:(Intercept)", ]
    expect_equal(ka$statistic, kb$statistic, tolerance = 1e-10)
    expect_equal(ka$p.value, kb$p.value, tolerance = 1e-10)
    expect_equal(ka$statistic, kc$statistic, tolerance = 1e-10)
  })

  it("effects = 'fixed' drops the variance rows", {
    fit <- .get_fit_for_methods()
    td <- tidy(fit, effects = "fixed")
    expect_true(all(td$component == "fixed"))
  })

  it("effects = 'ran_pars' returns only variance rows", {
    fit <- .get_fit_for_methods()
    td <- tidy(fit, effects = "ran_pars")
    expect_true(all(td$component == "variance"))
    expect_gt(nrow(td), 0L)
  })

  it("all four report_space values preserve the row count; log == internal (B11)", {
    fit <- .get_fit_for_methods()
    ni <- nrow(tidy(fit, report_space = "internal"))
    expect_equal(ni, nrow(tidy(fit, report_space = "natural")))
    expect_equal(ni, nrow(tidy(fit, report_space = "log10")))
    expect_equal(ni, nrow(tidy(fit, report_space = "log")))
    # B11: "log" is now requestable and coincides with "internal" for beta_k
    tl <- tidy(fit, report_space = "log", effects = "fixed")
    ti <- tidy(fit, report_space = "internal", effects = "fixed")
    expect_equal(tl$estimate, ti$estimate)
    expect_true(all(tl$estimate_scale == "log"))
  })

  it("emits an s shape row for a 2-parameter fit (8-col contract, natural)", {
    fit <- .get_gm_fit_for_methods()
    td <- tidy(fit, report_space = "natural")
    expect_named(td, c("term", "estimate", "std.error", "statistic",
                       "p.value", "component", "estimate_scale",
                       "term_display"))
    s_row <- td[td$term == "s", ]
    expect_equal(nrow(s_row), 1L)
    expect_equal(s_row$component, "shape")
    # natural-scale s == exp(log_s), with a real Wald SE (not NA)
    expect_equal(s_row$estimate,
                 exp(unname(fit$model$coefficients[["log_s"]])),
                 tolerance = 1e-8)
    expect_false(is.na(s_row$std.error))
  })

  it("internal space leaves s on the log scale (= log_s)", {
    fit <- .get_gm_fit_for_methods()
    td <- tidy(fit, report_space = "internal")
    s_row <- td[td$term == "s", ]
    expect_equal(s_row$estimate, unname(fit$model$coefficients[["log_s"]]),
                 tolerance = 1e-10)
    expect_equal(s_row$estimate_scale, "log")
  })

  it("a mazur fit emits NO s row", {
    fit <- .get_fit_for_methods()
    td <- tidy(fit)
    expect_false("s" %in% td$term)
  })
})

# ---------------------------------------------------------------------------
# M.6: glance
# ---------------------------------------------------------------------------

describe("glance", {
  it("returns one row with the canonical columns and backend string", {
    fit <- .get_fit_for_methods()
    g <- glance(fit)
    expect_s3_class(g, "tbl_df")
    expect_equal(nrow(g), 1L)
    expect_named(g, c("model_class", "backend", "equation", "family",
                      "nobs", "n_subjects", "n_random_effects",
                      "converged", "logLik", "AIC", "BIC"))
    expect_equal(g$model_class, "beezdiscounting_tmb")
    expect_equal(g$backend, "TMB_mixed")
    expect_equal(g$equation, "mazur")
    expect_equal(g$family, "sltb")
    expect_equal(g$nobs, fit$param_info$n_obs)
    expect_equal(g$n_subjects, fit$param_info$n_subjects)
    expect_equal(g$logLik, fit$loglik)
    expect_equal(g$AIC, fit$AIC)
    expect_equal(g$BIC, fit$BIC)
  })

  it("nobs / AIC / BIC / logLik are all numeric", {
    fit <- .get_fit_for_methods()
    g <- glance(fit)
    expect_type(g$nobs, "integer")
    expect_type(g$AIC, "double")
    expect_type(g$BIC, "double")
    expect_type(g$logLik, "double")
  })

  it("reports family = gaussian for a gaussian fit", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    d2 <- simulate_dd_ip(n_subjects = 25, family = "gaussian",
                                equation = "mazur", seed = 11)
    f2 <- fit_dd_tmb(d2, equation = "mazur", family = "gaussian",
                     random_effects = k ~ 1, verbose = 0)
    expect_equal(glance(f2)$family, "gaussian")
  })
})

# ---------------------------------------------------------------------------
# M.7: confint
# ---------------------------------------------------------------------------

describe("confint", {
  it("returns term/estimate/conf.low/conf.high/level with internal default", {
    fit <- .get_fit_for_methods()
    ci <- confint(fit)
    expect_s3_class(ci, "tbl_df")
    expect_named(ci, c("term", "estimate", "conf.low", "conf.high", "level"))
    expect_true(all(ci$conf.low <= ci$estimate))
    expect_true(all(ci$estimate <= ci$conf.high))
    expect_true(all(ci$level == 0.95))
  })

  it("Wald interval = estimate +/- z*se on the internal scale", {
    fit <- .get_fit_for_methods()
    ci <- confint(fit, parm = "k:(Intercept)")
    co <- fit$model$coefficients["beta_k"][1]
    se <- .dd_tmb_model_se(fit)["beta_k"]
    z <- stats::qnorm(0.975)
    expect_equal(ci$estimate, unname(co), tolerance = 1e-10)
    expect_equal(ci$conf.low, unname(co - z * se), tolerance = 1e-8)
    expect_equal(ci$conf.high, unname(co + z * se), tolerance = 1e-8)
  })

  it("report_space = 'natural' exponentiates beta_k bounds", {
    fit <- .get_fit_for_methods()
    ci_i <- confint(fit, parm = "k:(Intercept)", report_space = "internal")
    ci_n <- confint(fit, parm = "k:(Intercept)", report_space = "natural")
    expect_equal(ci_n$estimate, exp(ci_i$estimate), tolerance = 1e-8)
    expect_equal(ci_n$conf.low, exp(ci_i$conf.low), tolerance = 1e-8)
    expect_equal(ci_n$conf.high, exp(ci_i$conf.high), tolerance = 1e-8)
  })

  it("parm filters by display name OR raw name", {
    fit <- .get_fit_for_methods()
    ci_disp <- confint(fit, parm = "k:(Intercept)")
    ci_raw <- confint(fit, parm = "beta_k")
    expect_equal(nrow(ci_disp), 1L)
    expect_equal(nrow(ci_raw), 1L)
    expect_equal(ci_disp$estimate, ci_raw$estimate, tolerance = 1e-12)
  })

  it("honours a non-default confidence level", {
    fit <- .get_fit_for_methods()
    ci90 <- confint(fit, parm = "k:(Intercept)", level = 0.90)
    ci95 <- confint(fit, parm = "k:(Intercept)", level = 0.95)
    expect_true(ci90$conf.low > ci95$conf.low)
    expect_true(ci90$conf.high < ci95$conf.high)
    expect_true(all(ci90$level == 0.90))
  })

  it("includes an s row and back-transforms it under natural space", {
    fit <- .get_gm_fit_for_methods()
    ci_i <- confint(fit, parm = "s", report_space = "internal")
    ci_n <- confint(fit, parm = "s", report_space = "natural")
    expect_equal(nrow(ci_i), 1L)
    expect_equal(ci_n$estimate, exp(ci_i$estimate), tolerance = 1e-8)
    expect_equal(ci_n$conf.low, exp(ci_i$conf.low), tolerance = 1e-8)
    expect_equal(ci_n$conf.high, exp(ci_i$conf.high), tolerance = 1e-8)
    # finite interval => the s SE was populated (Codex 7)
    expect_true(is.finite(ci_i$conf.low) && is.finite(ci_i$conf.high))
  })

  it("filters the s row by raw name log_s too", {
    fit <- .get_gm_fit_for_methods()
    expect_equal(nrow(confint(fit, parm = "log_s")), 1L)
  })
})

# ---------------------------------------------------------------------------
# M.8: summary / print
# ---------------------------------------------------------------------------

describe("summary / print", {
  it("summary() returns the class with coefficients + variance components", {
    fit <- .get_fit_for_methods()
    s <- summary(fit)
    expect_s3_class(s, "summary.beezdiscounting_tmb")
    expect_equal(s$backend, "TMB_mixed")
    expect_equal(s$equation, "mazur")
    expect_equal(s$family, "sltb")
    expect_s3_class(s$coefficients, "tbl_df")
    expect_true(all(c("term", "estimate", "std.error", "statistic",
                      "p.value") %in% names(s$coefficients)))
    expect_true(!is.null(s$variance_components))
    expect_equal(s$n_subjects, fit$param_info$n_subjects)
  })

  it("summary() honours report_space for the fixed-effect estimates", {
    fit <- .get_fit_for_methods()
    s_nat <- summary(fit, report_space = "natural")
    krow <- s_nat$coefficients[s_nat$coefficients$term == "k:(Intercept)", ]
    expect_equal(krow$estimate,
                 exp(unname(fit$model$coefficients["beta_k"][1])),
                 tolerance = 1e-8)
  })

  it("summary(report_space = 'log') == 'internal' for beta_k (B11)", {
    fit <- .get_fit_for_methods()
    sl <- summary(fit, report_space = "log")
    si <- summary(fit, report_space = "internal")
    expect_equal(sl$coefficients$estimate, si$coefficients$estimate)
    krow <- sl$coefficients[sl$coefficients$term == "k:(Intercept)", ]
    expect_equal(krow$estimate, unname(fit$model$coefficients["beta_k"][1]),
                 tolerance = 1e-8)
  })

  it("print() and print.summary() run without error and return invisibly", {
    fit <- .get_fit_for_methods()
    expect_invisible(print(fit))
    expect_output(print(fit), "TMB Mixed-Effects Discounting Model")
    s <- summary(fit)
    expect_invisible(print(s))
    expect_output(print(s), "TMB_mixed")
  })

  it("print.summary header reflects the report space / scale (R4)", {
    fit <- .get_fit_for_methods()
    # natural -> "(k)"
    expect_output(print(summary(fit, report_space = "natural")),
                  "Fixed Effects \\(k\\)")
    # log10 -> "(log10 k)"
    expect_output(print(summary(fit, report_space = "log10")),
                  "Fixed Effects \\(log10 k\\)")
    # internal -> "(log k)"
    expect_output(print(summary(fit, report_space = "internal")),
                  "Fixed Effects \\(log k\\)")
  })

  it("summary notes flag non-convergence / missing SEs", {
    fit <- .get_fit_for_methods()
    s <- summary(fit)
    expect_type(s$notes, "character")
  })

  it("summary() carries the fitting call, not the optimizer status string", {
    fit <- .get_fit_for_methods()
    s <- summary(fit)
    expect_identical(s$call, fit$call)
    # it must NOT be the optimizer message (R5)
    expect_false(identical(s$call, fit$opt$message))
    expect_true(is.call(s$call) || is.null(s$call))
  })

  it("summary() shows s in the coefficient table for a 2-parameter fit", {
    fit <- .get_gm_fit_for_methods()
    s <- summary(fit, report_space = "natural")
    srow <- s$coefficients[s$coefficients$term == "s", ]
    expect_equal(nrow(srow), 1L)
    expect_equal(srow$component, "shape")
    expect_equal(srow$estimate,
                 exp(unname(fit$model$coefficients[["log_s"]])),
                 tolerance = 1e-8)
  })

  it("summary() of a mazur fit has no s row", {
    fit <- .get_fit_for_methods()
    s <- summary(fit)
    expect_false("s" %in% s$coefficients$term)
  })

  it("summary notes flag the shape parameter for a 2-parameter fit", {
    fit <- .get_gm_fit_for_methods()
    s <- summary(fit)
    expect_true(any(grepl("shape parameter", s$notes)))
    expect_output(print(s), "shape parameter")
  })
})

describe("non-default column names (B1: canonical-names contract)", {
  it("predict/fitted/residuals/augment work when y_var/x_var/id_var are remapped", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 18, seed = 5)
    remap <- sim
    names(remap) <- c("subj", "delay", "indiff")  # id, x, y -> remapped names

    fit_canon <- fit_dd_tmb(sim, equation = "mazur", family = "sltb", verbose = 0)
    fit_remap <- fit_dd_tmb(remap, y_var = "indiff", x_var = "delay",
                            id_var = "subj", equation = "mazur",
                            family = "sltb", verbose = 0)

    # Regression: these previously errored because param_info kept the user's
    # original names while fit$data is canonical id/x/y.
    expect_no_error(p <- predict(fit_remap, type = "response"))
    expect_true(".fitted" %in% names(p))
    expect_equal(length(fitted(fit_remap)), fit_remap$param_info$n_obs)
    expect_equal(length(residuals(fit_remap)), fit_remap$param_info$n_obs)
    aug <- augment(fit_remap)
    expect_true(all(c(".fitted", ".resid", ".std_resid") %in% names(aug)))

    # Identical underlying data (+ deterministic optimizer) -> identical fitted.
    expect_equal(unname(fitted(fit_remap)), unname(fitted(fit_canon)),
                 tolerance = 1e-6)
  })

  it("predict() errors cleanly when newdata omits the canonical delay column", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    fit <- .get_fit_for_methods()
    bad <- data.frame(id = "s1", y = 0.5)   # canonical 'x' (delay) absent
    expect_error(
      predict(fit, newdata = bad, type = "response"),
      "delay column|canonical"
    )
  })
})

describe("confint se_available gate (B2)", {
  it("returns NA intervals (keeps estimates) and warns when se_available is FALSE", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    fit <- .get_fit_for_methods()
    fit$se_available <- FALSE
    expect_warning(ci <- confint(fit), "unreliable")
    expect_true(all(is.na(ci$conf.low)) && all(is.na(ci$conf.high)))
    expect_false(any(is.na(ci$estimate)))   # point estimates preserved
  })

  it("finite CIs on a normal (PD-Hessian) fit (no over-broadening)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    fit <- .get_fit_for_methods()
    if (isTRUE(fit$se_available)) {
      ci <- confint(fit)
      expect_true(any(is.finite(ci$conf.low)))
    } else {
      succeed("fixture fit had unavailable SEs; gate-on path covered above")
    }
  })
})

describe("2-RE (k + phi) S3 surface", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  make_fit <- function(cov = "pdSymm", seed = 21) {
    key <- paste0("re2_", cov, "_", seed)
    if (!exists(key, envir = .ddm_cache)) {
      sim <- simulate_dd_ip(n_subjects = 30, sigma_u = 0.5, sigma_phi = 0.4,
                            rho_kphi = 0.3, phi = 10, seed = seed)
      .ddm_cache[[key]] <- fit_dd_tmb(sim, random_effects = k + phi ~ 1,
                                      covariance_structure = cov, verbose = 0)
    }
    .ddm_cache[[key]]
  }

  it("VarCorr returns a 2x2 structure with SDs and a correlation (pdSymm)", {
    vc <- VarCorr(make_fit("pdSymm"))
    expect_true(is.data.frame(vc) || is.matrix(vc))
    # two variance rows + one correlation entry expected
    expect_gte(NROW(vc), 2L)
  })

  it("VarCorr correlation is a structural 0 for pdDiag", {
    vc <- VarCorr(make_fit("pdDiag"))
    # pdDiag fixes cor_re out -> rho = tanh(0) = 0, so Sigma is diagonal and
    # VarCorr reports Corr = c(NA, 0) (first row's Corr is NA by convention).
    expect_true(is.na(vc$Corr[1]))
    expect_equal(vc$Corr[2], 0, tolerance = 1e-10)
  })

  it("ranef returns per-subject k and phi for a 2-RE fit", {
    re <- ranef(make_fit("pdSymm"))
    expect_true(all(c("id", "k", "phi") %in% names(re)))
    expect_gt(length(unique(round(re$phi, 6))), 1L)
  })

  it("ranef for a 1-RE fit is unchanged (no phi column)", {
    sim <- simulate_dd_ip(n_subjects = 12, seed = 5)
    re <- ranef(fit_dd_tmb(sim, verbose = 0))
    expect_false("phi" %in% names(re))
    expect_true(all(c("id", "u_i", "k") %in% names(re)))
  })

  it("predict/fitted/residuals/augment run on a 2-RE fit", {
    fit <- make_fit("pdSymm")
    # predict() returns a tibble with one row per observation plus .fitted.
    pred <- predict(fit)
    expect_equal(nrow(pred), nrow(fit$data))
    expect_true(all(is.finite(pred$.fitted)))
    expect_true(all(is.finite(fitted(fit))))
    expect_true(all(is.finite(residuals(fit, type = "pearson"))))
    aug <- augment(fit)
    expect_true(all(c(".fitted", ".resid", ".std_resid") %in% names(aug)))
    expect_true(all(is.finite(aug$.std_resid)))
  })
})
