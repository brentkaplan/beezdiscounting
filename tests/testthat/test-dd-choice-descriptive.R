skip_on_cran(); skip_if_not_installed("TMB")

describe("fit_dd_choice(mode='descriptive') — pipeline + map", {
  it("fits, converges, and is classed beezdiscounting_choice / descriptive", {
    dat <- .choice_desc_fixture(n_subjects = 60, seed = 12)
    fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
    expect_s3_class(fit, "beezdiscounting_choice")
    expect_equal(fit$param_info$mode, "descriptive")
    expect_equal(fit$param_info$n_random_effects, 2L)
    expect_true(fit$converged)
  })

  it("free parameters match the descriptive mode (no structural params free)", {
    dat <- .choice_desc_fixture(n_subjects = 50, seed = 13)
    fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
    free <- names(fit$opt$par)
    expect_true(all(c("theta", "log_sd_re", "cor_re") %in% free))
    expect_false(any(c("beta_k", "log_gamma", "beta0", "log_sigma_u") %in% free))
    expect_equal(length(fit$opt$par), 5L)   # theta(2)+log_sd_re(2)+cor_re(1); b integrated
  })

  it("random_slopes = FALSE => pooled (no RE params free), q = 0", {
    dat <- .choice_desc_fixture(n_subjects = 50, seed = 14)
    fit <- fit_dd_choice(dat, mode = "descriptive", random_slopes = FALSE, verbose = 0)
    free <- names(fit$opt$par)
    expect_true("theta" %in% free)
    expect_false(any(c("log_sd_re", "cor_re") %in% free))
    expect_equal(fit$param_info$n_random_effects, 0L)
  })

  it("predictors override extends the fixed design (Z gains a column)", {
    dat <- .choice_desc_fixture(n_subjects = 50, seed = 15)
    dat$grp <- factor(rep(c("x", "y"), length.out = nrow(dat)))
    fit <- fit_dd_choice(
      dat, mode = "descriptive",
      predictors = ~ 0 + log(ll_amount / ss_amount) + log(delay + 1) + grp,
      verbose = 0)
    expect_gte(sum(names(fit$opt$par) == "theta"), 3L)  # 2 Young + grp dummies
    expect_equal(fit$param_info$n_random_effects, 2L)   # Zre unchanged
  })

  it("structural mode still fits + no descriptive param leaks free (regression)", {
    dat <- .choice_fit_fixture(n_subjects = 30, seed = 9)
    fit <- fit_dd_choice(dat, mode = "structural", equation = "mazur", verbose = 0)
    expect_true(fit$converged)
    expect_equal(fit$param_info$mode, "structural")
    free <- names(fit$opt$par)
    expect_true(all(c("beta_k", "log_sigma_u", "log_gamma") %in% free))
    expect_false(any(c("theta", "log_sd_re", "cor_re", "beta0") %in% free))
  })

  it("warns that factors/continuous_covariates are ignored for descriptive", {
    dat <- .choice_desc_fixture(n_subjects = 30, seed = 16)
    dat$grp <- factor(rep(c("x", "y"), length.out = nrow(dat)))
    expect_warning(
      fit_dd_choice(dat, mode = "descriptive", factors = "grp", verbose = 0),
      regexp = "ignored|descriptive|predictors"
    )
  })
})
