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
    d2 <- .simulate_dd_ip_mixed(n_subjects = 25, family = "gaussian",
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
    d3 <- .simulate_dd_ip_mixed(n_subjects = 30, family = "sltb",
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
    d2 <- .simulate_dd_ip_mixed(n_subjects = 25, family = "gaussian",
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
