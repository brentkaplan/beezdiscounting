# P1 contract guards from the Phase-1 audit, tranche 2 (B4-B8).

describe("between-subject contract warning (B4)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("warns when a declared covariate varies within a subject", {
    d <- simulate_dd_ip(n_subjects = 12, seed = 4)
    d$age <- as.numeric(d$x)   # varies within subject (= the delay)
    expect_warning(
      suppressMessages(fit_dd_tmb(d, continuous_covariates = "age", verbose = 0)),
      "not constant within"
    )
  })

  it("does not warn for a genuinely between-subject covariate", {
    d <- simulate_dd_ip(n_subjects = 12, seed = 5)
    ids <- levels(d$id)
    set.seed(1)
    d$age <- stats::setNames(stats::rnorm(length(ids), 30, 5), ids)[as.character(d$id)]
    w <- testthat::capture_warnings(
      suppressMessages(fit_dd_tmb(d, continuous_covariates = "age", verbose = 0))
    )
    expect_false(any(grepl("not constant within", w)))
  })

  it("warns (and counts subjects) for a FACTOR that varies within a subject", {
    d <- simulate_dd_ip(n_subjects = 10, seed = 6)
    d$grp <- factor(rep(c("p", "q"), length.out = nrow(d)))  # varies within subject
    expect_warning(
      suppressMessages(fit_dd_tmb(d, factors = "grp", verbose = 0)),
      "not constant within"
    )
  })
})

describe("declared-factor coercion + covariate type (B5)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("coerces a numeric column passed as factors= to dummy columns, not a slope", {
    d <- simulate_dd_ip(n_subjects = 12, n_conditions = 2, delta_k = c(0, 0.5),
                        seed = 2)
    d$grp <- as.integer(d$condition)   # numeric 1/2
    fit <- fit_dd_tmb(d, factors = "grp", verbose = 0)
    # A factor -> dummy column grp2; a continuous slope would be a single 'grp'.
    expect_true("grp2" %in% colnames(fit$formula_details$X))
    expect_false("grp" %in% colnames(fit$formula_details$X))
  })

  it("rejects a non-finite continuous covariate", {
    d <- simulate_dd_ip(n_subjects = 10, seed = 3)
    d$age <- 30
    d$age[1] <- Inf
    expect_error(fit_dd_tmb(d, continuous_covariates = "age", verbose = 0), "finite")
  })

  it("rejects a non-numeric continuous covariate", {
    d <- simulate_dd_ip(n_subjects = 10, seed = 3)
    d$age <- rep(c("low", "high"), length.out = nrow(d))   # character covariate
    expect_error(fit_dd_tmb(d, continuous_covariates = "age", verbose = 0),
                 "numeric")
  })

  it("rejects a predictor that names the id / delay / response column", {
    d <- simulate_dd_ip(n_subjects = 10, seed = 1)
    expect_error(fit_dd_tmb(d, factors = "x", verbose = 0), "cannot be the id")
  })
})

describe("reserved canonical names (B6)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("rejects an extra column named id/x/y under a remapped role var", {
    d <- simulate_dd_ip(n_subjects = 10, seed = 1)
    d$pid <- d$id
    d$id  <- factor(rep(c("g1", "g2"), length.out = nrow(d)))  # a factor literally 'id'
    expect_error(
      fit_dd_tmb(d, id_var = "pid", factors = "id", verbose = 0),
      "reserved"
    )
  })
})

describe("generalized blow-up predicate (B7)", {
  it(".dd_logk_blowup uses the full eta = X beta_k, not only the intercept", {
    X <- cbind(1, c(0, 1, 0, 1))
    # sane intercept + sane non-reference coefficient -> not a blow-up
    expect_false(.dd_logk_blowup(list(par = c(beta_k = log(0.01), beta_k = 0.5)), X))
    # sane intercept but a blown-up NON-reference coefficient -> blow-up (the
    # case the old beta_k[1]-only guard missed)
    expect_true(.dd_logk_blowup(list(par = c(beta_k = log(0.01), beta_k = 30)), X))
    # non-finite beta and beta/X shape mismatch -> blow-up
    expect_true(.dd_logk_blowup(list(par = c(beta_k = NA_real_, beta_k = 1)), X))
    expect_true(.dd_logk_blowup(list(par = c(beta_k = 0)), X))
    # non-finite eta (a non-finite design entry) -> blow-up
    expect_true(.dd_logk_blowup(list(par = c(beta_k = 1, beta_k = 1)),
                                cbind(1, c(0, Inf))))
    # intercept-only: max|X beta| == |beta[1]|, so behavior is unchanged
    Xi <- cbind(rep(1, 4))
    expect_false(.dd_logk_blowup(list(par = c(beta_k = 5)), Xi))
    expect_true(.dd_logk_blowup(list(par = c(beta_k = 25)), Xi))
  })
})

describe("predict newdata column guard (B8)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("errors cleanly when newdata omits a factor/covariate column", {
    d <- simulate_dd_ip(n_subjects = 12, n_conditions = 2, delta_k = c(0, 0.5),
                        seed = 2)
    fit <- fit_dd_tmb(d, factors = "condition", verbose = 0)
    expect_error(
      predict(fit, newdata = data.frame(x = c(7, 365))),   # omits 'condition'
      "missing predictor"
    )
  })

  it("works when the factor column is supplied (population level)", {
    d <- simulate_dd_ip(n_subjects = 12, n_conditions = 2, delta_k = c(0, 0.5),
                        seed = 2)
    fit <- fit_dd_tmb(d, factors = "condition", verbose = 0)
    nd  <- data.frame(x = c(7, 365),
                      condition = factor("C1", levels = c("C1", "C2")))
    expect_no_error(predict(fit, newdata = nd, level = "population"))
  })
})
