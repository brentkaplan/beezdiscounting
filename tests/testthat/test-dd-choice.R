describe(".dd_choice_build_map()", {
  starts <- list(
    beta_k = 0, log_sigma_u = 0, log_gamma = 0, beta0 = 0,
    u = matrix(0, 2L, 1L),
    theta = 0, log_sd_re = rep(0, 2L), cor_re = 0, b = matrix(0, 2L, 2L)
  )
  it("always fixes the descriptive blocks (theta/log_sd_re/cor_re/b)", {
    for (intercept in c(TRUE, FALSE)) {
      m <- .dd_choice_build_map(intercept = intercept, starts = starts)
      expect_true(all(c("theta", "log_sd_re", "cor_re", "b") %in% names(m)))
      expect_true(all(is.na(m$theta)))
      expect_equal(length(m$log_sd_re), length(starts$log_sd_re))
      expect_true(all(is.na(m$log_sd_re)))
      expect_true(all(is.na(m$cor_re)))
      expect_equal(length(m$b), length(starts$b))
      expect_true(all(is.na(m$b)))
    }
  })
  it("maps beta0 when intercept is FALSE and frees it when TRUE", {
    m_off <- .dd_choice_build_map(intercept = FALSE, starts = starts)
    expect_true("beta0" %in% names(m_off)); expect_true(all(is.na(m_off$beta0)))
    expect_false("beta0" %in% names(.dd_choice_build_map(intercept = TRUE,
                                                         starts = starts)))
  })
})

describe("fit_dd_choice() structural", {
  it("fits, returns the beezdiscounting_choice contract (intercept off)", {
    skip_on_cran(); skip_if_not_installed("TMB")
    dat <- .choice_fit_fixture(seed = 11)
    fit <- fit_dd_choice(dat, mode = "structural", equation = "mazur", verbose = 0)
    expect_s3_class(fit, "beezdiscounting_choice")
    expect_true(all(c("call", "opt", "model", "sdr", "param_info",
                      "formula_details", "subject_pars", "loglik", "AIC", "BIC",
                      "converged", "se_available", "data") %in% names(fit)))
    co <- fit$model$coefficients
    expect_equal(names(fit$opt$par)[1], "beta_k")
    expect_false("beta0" %in% names(fit$opt$par))
    expect_true("log_gamma" %in% names(co))
    expect_equal(fit$param_info$mode, "structural")
    expect_equal(fit$param_info$equation, "mazur")
    expect_false(fit$param_info$intercept)
    expect_named(fit$subject_pars, c("id", "u_i", "k"))
    expect_equal(nrow(fit$subject_pars), 40L)
  })

  it("includes a free beta0 when intercept = TRUE (df + 1)", {
    skip_on_cran(); skip_if_not_installed("TMB")
    dat <- .choice_fit_fixture(beta0 = 0.5, seed = 12)
    off <- fit_dd_choice(dat, mode = "structural", intercept = FALSE, verbose = 0)
    on  <- fit_dd_choice(dat, mode = "structural", intercept = TRUE, verbose = 0)
    expect_false("beta0" %in% names(off$opt$par))
    expect_true("beta0" %in% names(on$opt$par))
    expect_equal(length(on$opt$par), length(off$opt$par) + 1L)
  })

  it("rejects descriptive mode (Plan B) cleanly", {
    skip_on_cran(); skip_if_not_installed("TMB")
    dat <- .choice_fit_fixture(seed = 13)
    expect_error(fit_dd_choice(dat, mode = "descriptive", verbose = 0),
                 "descriptive|not yet|Plan B|implemented")
  })
})
