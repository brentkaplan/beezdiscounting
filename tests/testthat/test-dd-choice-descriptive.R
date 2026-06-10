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

# Well-conditioned descriptive design: 6 magnitude pairs (ll/ss ratios 1.25-4)
# crossed with 5 delays -> strong variation in BOTH log(ll/ss) and log(delay+1),
# so the 2-RE (co)variances are identifiable. (The MCQ27 stimulus set has ll/ss
# ratios near 1 and only weakly identifies the magnitude-slope variance, so it is
# used below only for sign-level / plumbing checks, not (co)variance recovery.)
.choice_wellcond_design <- function() {
  mp <- list(c(40, 50), c(30, 60), c(25, 75), c(50, 80), c(20, 60), c(15, 60))
  del <- c(3, 14, 45, 120, 300)
  g <- expand.grid(i = seq_along(mp), d = del)
  list(ss = vapply(g$i, function(j) mp[[j]][1], numeric(1)),
       ll = vapply(g$i, function(j) mp[[j]][2], numeric(1)),
       delay = g$d)
}

describe("fit_dd_choice(mode='descriptive') — recovery", {
  it("recovers theta and Sigma from known-truth data (centered seed 308)", {
    # Seed 308 chosen from a 12-seed sweep (301-312) as CENTERED on all five
    # quantities (not a barely-passing edge): recovered theta 1.231/-0.622,
    # sd 0.401/0.306, rho -0.244 vs truth theta 1.2/-0.6, sd 0.5/0.3, rho -0.3;
    # all swept seeds converged with a PD Hessian. The magnitude-slope SD is the
    # weakest-identified cell, hence its looser absolute tolerance.
    dz <- .choice_wellcond_design()
    dat <- simulate_dd_choice(
      n_subjects = 150, mode = "descriptive",
      ss_amount = dz$ss, ll_amount = dz$ll, delay = dz$delay,
      theta = c(1.2, -0.6), re_sd = c(0.5, 0.3), re_cor = -0.3, seed = 308)
    fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
    expect_true(fit$converged)
    th <- unname(fit$model$coefficients[names(fit$model$coefficients) == "theta"])
    expect_equal(th[1],  1.2, tolerance = 0.25)
    expect_equal(th[2], -0.6, tolerance = 0.20)
    sds <- sqrt(diag(fit$Sigma))
    expect_equal(sds[1], 0.5, tolerance = 0.25)   # magnitude-slope SD (loosest)
    expect_equal(sds[2], 0.3, tolerance = 0.12)   # delay-slope SD (well-identified)
    rho <- fit$Sigma[1, 2] / prod(sds)
    expect_equal(rho, -0.3, tolerance = 0.30)
  })

  it("recovers Young's sign pattern (mag > 0, delay < 0)", {
    dz <- .choice_wellcond_design()
    dat <- simulate_dd_choice(
      n_subjects = 120, mode = "descriptive",
      ss_amount = dz$ss, ll_amount = dz$ll, delay = dz$delay,
      theta = c(1.4, -0.6), re_sd = c(0.4, 0.25), re_cor = -0.2, seed = 309)
    fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
    th <- unname(fit$model$coefficients[names(fit$model$coefficients) == "theta"])
    expect_gt(th[1], 0); expect_lt(th[2], 0)
  })
})

describe("fit_dd_choice(mode='descriptive') — MCQ27 integration", {
  it("mcq27_to_choice() output fits as a descriptive model (plumbing smoke)", {
    # generate_data_mcq() responses are random coin flips (no amount/delay
    # signal) => this is a PLUMBING smoke test of the adapter -> fit path, NOT a
    # recovery claim. Assert only fit/class/finite.
    set.seed(202)
    resp <- generate_data_mcq(n_ids = 60, n_items = 27)
    choice_frame <- mcq27_to_choice(resp)
    fit <- fit_dd_choice(choice_frame, mode = "descriptive", verbose = 0)
    expect_s3_class(fit, "beezdiscounting_choice")
    th <- unname(fit$model$coefficients[names(fit$model$coefficients) == "theta"])
    expect_true(all(is.finite(th)))
  })

  it("recovers Young's sign pattern on the REAL MCQ27 stimulus design", {
    # Real 27-item Kirby (ss, ll, delay) design with choices from a KNOWN
    # descriptive truth. Theta SIGNS are robustly identified on this design even
    # though the magnitude-slope (co)variances are not (ll/ss ratios near 1).
    lt <- get_lookup_table()
    dat <- simulate_dd_choice(
      n_subjects = 120, mode = "descriptive",
      ss_amount = lt$ss_amount, ll_amount = lt$ll_amount, delay = lt$delay,
      theta = c(1.5, -0.5), re_sd = c(0.4, 0.25), re_cor = -0.2, seed = 203)
    fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
    th <- unname(fit$model$coefficients[names(fit$model$coefficients) == "theta"])
    expect_gt(th[1], 0); expect_lt(th[2], 0)
  })
})
