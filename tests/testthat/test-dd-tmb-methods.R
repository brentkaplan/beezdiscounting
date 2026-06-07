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
