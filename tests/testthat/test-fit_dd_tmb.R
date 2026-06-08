describe(".dd_tmb_prepare_data()", {
  it("0-indexes subject_id and aligns vectors with subject_levels", {
    dat <- data.frame(
      id = factor(c("b", "b", "a", "a"), levels = c("b", "a")),
      x = c(7, 30, 7, 30),
      y = c(0.9, 0.5, 0.8, 0.4)
    )
    prep <- .dd_tmb_prepare_data(dat, y_var = "y", x_var = "x", id_var = "id")
    expect_equal(prep$subject_levels, c("a", "b"))
    expect_equal(prep$n_subjects, 2L)
    expect_equal(prep$n_obs, 4L)
    expect_equal(min(prep$subject_id), 0L)
    expect_equal(max(prep$subject_id), 1L)
    # subject_id maps each row to the sorted-level index
    expect_equal(prep$subject_id, c(1L, 1L, 0L, 0L))
    expect_type(prep$subject_id, "integer")
  })

  it("drops rows with NA in y/x/id and keeps x,y,subject_id row-coherent", {
    dat <- data.frame(
      id = c("a", "a", "b", "b"),
      x = c(7, 30, 7, NA),
      y = c(0.9, NA, 0.8, 0.4)
    )
    prep <- .dd_tmb_prepare_data(dat, y_var = "y", x_var = "x", id_var = "id")
    expect_equal(prep$n_obs, 2L)
    expect_equal(prep$y, c(0.9, 0.8))
    expect_equal(prep$x, c(7, 7))
    expect_equal(nrow(prep$data), 2L)
    # surviving rows keep their parallel arrays aligned
    expect_equal(length(prep$x), length(prep$subject_id))
  })

  it("renames caller columns to canonical id/x/y in $data", {
    dat <- data.frame(pid = "p1", delay = 7, ip = 0.5)
    prep <- .dd_tmb_prepare_data(dat, y_var = "ip", x_var = "delay", id_var = "pid")
    expect_true(all(c("id", "x", "y") %in% names(prep$data)))
    expect_equal(names(prep$data)[1:3], c("id", "x", "y"))
  })

  it("retains extra cols and complete-cases id/x/y + extras ONCE (row-coherent)", {
    dat <- data.frame(
      id = c("a", "a", "b", "b"),
      x = c(7, 30, 7, 30),
      y = c(0.9, 0.5, 0.8, 0.4),
      grp = factor(c("ctrl", "ctrl", NA, "trt")),  # one NA in an extra col
      age = c(20, 20, 30, 30),
      stringsAsFactors = FALSE
    )
    prep <- .dd_tmb_prepare_data(dat, "y", "x", "id",
      extra_cols = c("grp", "age"))
    # the row with NA grp is dropped because grp is a modeling column
    expect_equal(prep$n_obs, 3L)
    # $data keeps id/x/y first then the retained extras, all row-aligned
    expect_equal(names(prep$data), c("id", "x", "y", "grp", "age"))
    expect_equal(prep$data$y, c(0.9, 0.5, 0.4))
    expect_equal(as.character(prep$data$grp), c("ctrl", "ctrl", "trt"))
    # parallel arrays match $data exactly (one filtered row order)
    expect_equal(prep$y, prep$data$y)
    expect_equal(prep$x, prep$data$x)
    expect_equal(length(prep$subject_id), nrow(prep$data))
    # unused factor levels dropped on the surviving frame
    expect_setequal(levels(prep$data$grp), c("ctrl", "trt"))
  })
})

describe(".dd_tmb_build_design()", {
  it("returns an intercept-only X (+ rhs formula + contrasts) with no factors", {
    dat <- data.frame(id = "a", x = 7, y = 0.5)
    d <- .dd_tmb_build_design(dat)
    expect_equal(colnames(d$X), "(Intercept)")
    expect_equal(ncol(d$X), 1L)
    # rhs is a one-sided formula from build_fixed_rhs(), not a string
    expect_s3_class(d$rhs, "formula")
    expect_equal(rlang::f_text(d$rhs), "1")
    # contrasts attr is carried through (NULL when no factors)
    expect_true("contrasts" %in% names(d))
  })

  it("expands a single between-subject factor into contrast columns", {
    dat <- data.frame(
      id = rep(c("a", "b", "c", "d"), each = 1),
      x = 7,
      y = 0.5,
      grp = factor(c("ctrl", "ctrl", "trt", "trt"))
    )
    d <- .dd_tmb_build_design(dat, factors = "grp")
    expect_equal(ncol(d$X), 2L)
    expect_true("grptrt" %in% colnames(d$X))
    expect_match(rlang::f_text(d$rhs), "grp")
    # the per-factor contrast scheme is stored for the emmeans ref grid (E.1)
    expect_true("grp" %in% names(d$contrasts))
  })

  it("adds an interaction term when factor_interaction = TRUE", {
    dat <- data.frame(
      id = letters[1:4], x = 7, y = 0.5,
      a = factor(c("x", "x", "z", "z")),
      b = factor(c("p", "q", "p", "q"))
    )
    d <- .dd_tmb_build_design(dat, factors = c("a", "b"),
                              factor_interaction = TRUE)
    expect_true(any(grepl(":", colnames(d$X))))
  })

  it("adds continuous covariate main-effect columns", {
    dat <- data.frame(id = letters[1:3], x = 7, y = 0.5, age = c(20, 30, 40))
    d <- .dd_tmb_build_design(dat, continuous_covariates = "age")
    expect_true("age" %in% colnames(d$X))
  })
})

describe(".dd_tmb_build_tmb_data()", {
  it("assembles the TMB data list with the contract names and enum ints", {
    prep <- list(
      y = c(0.9, 0.5), x = c(7, 30), subject_id = c(0L, 0L),
      subject_levels = "a", n_subjects = 1L, n_obs = 2L,
      data = data.frame(id = "a", x = c(7, 30), y = c(0.9, 0.5))
    )
    design <- list(X = matrix(1, nrow = 2, ncol = 1,
                              dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    td <- .dd_tmb_build_tmb_data(prep, design, equation = "mazur",
                                 family = "sltb")
    expect_equal(td$model, "MixedDiscounting")
    expect_equal(td$eqn_type, 0L)   # mazur
    expect_equal(td$family, 0L)     # sltb
    expect_equal(td$n_obs, 2L)
    expect_equal(td$n_subjects, 1L)
    expect_equal(td$subject_id, c(0L, 0L))
    expect_true(is.matrix(td$X))
  })

  it("maps exponential -> 1 and gaussian -> 1", {
    prep <- list(y = 0.5, x = 7, subject_id = 0L, subject_levels = "a",
                 n_subjects = 1L, n_obs = 1L,
                 data = data.frame(id = "a", x = 7, y = 0.5))
    design <- list(X = matrix(1, 1, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    td <- .dd_tmb_build_tmb_data(prep, design, equation = "exponential",
                                 family = "gaussian")
    expect_equal(td$eqn_type, 1L)
    expect_equal(td$family, 1L)
  })

  it("maps green-myerson -> 2 and rachlin -> 3", {
    prep <- list(y = 0.5, x = 7, subject_id = 0L, subject_levels = "a",
                 n_subjects = 1L, n_obs = 1L,
                 data = data.frame(id = "a", x = 7, y = 0.5))
    design <- list(X = matrix(1, 1, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    gm <- .dd_tmb_build_tmb_data(prep, design, equation = "green-myerson",
                                 family = "sltb")
    ra <- .dd_tmb_build_tmb_data(prep, design, equation = "rachlin",
                                 family = "sltb")
    expect_equal(gm$eqn_type, 2L)
    expect_equal(ra$eqn_type, 3L)
  })
})

describe(".dd_tmb_default_starts()", {
  it("derives log_k intercept from median y at min delay (mazur)", {
    # median y at min delay = 0.8 -> mu=0.8 -> k = (1/mu - 1)/x_min
    prep <- list(
      y = c(0.8, 0.8, 0.3, 0.3), x = c(7, 7, 365, 365),
      subject_id = c(0L, 1L, 0L, 1L), subject_levels = c("a", "b"),
      n_subjects = 2L, n_obs = 4L,
      data = data.frame(id = c("a", "b", "a", "b"),
                        x = c(7, 7, 365, 365), y = c(0.8, 0.8, 0.3, 0.3))
    )
    design <- list(X = matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb")
    expect_length(st$beta_k, 1L)
    k_implied <- (1 / 0.8 - 1) / 7
    expect_equal(st$beta_k[1], log(k_implied), tolerance = 1e-8)
    expect_equal(st$log_sigma_u, log(0.5))
    expect_equal(st$log_aux, log(8))          # sltb
    expect_equal(dim(st$u), c(2L, 1L))
  })

  it("uses log(0.1) for log_aux under gaussian", {
    prep <- list(y = c(0.8, 0.3), x = c(7, 365), subject_id = c(0L, 0L),
                 subject_levels = "a", n_subjects = 1L, n_obs = 2L,
                 data = data.frame(id = "a", x = c(7, 365), y = c(0.8, 0.3)))
    design <- list(X = matrix(1, 2, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "gaussian")
    expect_equal(st$log_aux, log(0.1))
  })

  it("derives log_k intercept from exponential inversion (equation = 'exponential')", {
    # median y at min delay (x=7) = 0.8 -> mu=0.8 -> k = -log(mu)/x_min
    prep <- list(
      y = c(0.8, 0.8, 0.3, 0.3), x = c(7, 7, 365, 365),
      subject_id = c(0L, 1L, 0L, 1L), subject_levels = c("a", "b"),
      n_subjects = 2L, n_obs = 4L,
      data = data.frame(id = c("a", "b", "a", "b"),
                        x = c(7, 7, 365, 365), y = c(0.8, 0.8, 0.3, 0.3))
    )
    design <- list(X = matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb",
                                 equation = "exponential")
    k_implied <- -log(0.8) / 7
    expect_equal(st$beta_k[1], log(k_implied), tolerance = 1e-8)
  })

  it("falls back to log(0.01) intercept when x_min <= 0", {
    # x_min = 0: inversion would divide by zero; must use safe default k0 = 0.01
    prep <- list(
      y = c(0.9, 0.8, 0.3, 0.3), x = c(0, 0, 365, 365),
      subject_id = c(0L, 1L, 0L, 1L), subject_levels = c("a", "b"),
      n_subjects = 2L, n_obs = 4L,
      data = data.frame(id = c("a", "b", "a", "b"),
                        x = c(0, 0, 365, 365), y = c(0.9, 0.8, 0.3, 0.3))
    )
    design <- list(X = matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st_mazur <- .dd_tmb_default_starts(prep, design, family = "sltb",
                                       equation = "mazur")
    st_exp <- .dd_tmb_default_starts(prep, design, family = "sltb",
                                     equation = "exponential")
    expect_equal(st_mazur$beta_k[1], log(0.01), tolerance = 1e-8)
    expect_equal(st_exp$beta_k[1],   log(0.01), tolerance = 1e-8)
  })

  it("zero-pads beta_k for multi-column designs", {
    prep <- list(y = c(0.8, 0.5), x = c(7, 7), subject_id = c(0L, 1L),
                 subject_levels = c("a", "b"), n_subjects = 2L, n_obs = 2L,
                 data = data.frame(id = c("a", "b"), x = 7, y = c(0.8, 0.5)))
    design <- list(
      X = matrix(c(1, 1, 0, 1), 2, 2,
                 dimnames = list(NULL, c("(Intercept)", "grptrt"))),
      rhs = "~ grp")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb")
    expect_length(st$beta_k, 2L)
    expect_equal(st$beta_k[2], 0)
  })

  it("always includes log_s = 0 (s = 1 start; map-fixed for 1-parameter eqns)", {
    prep <- list(y = c(0.8, 0.3), x = c(7, 365), subject_id = c(0L, 0L),
                 subject_levels = "a", n_subjects = 1L, n_obs = 2L,
                 data = data.frame(id = "a", x = c(7, 365), y = c(0.8, 0.3)))
    design <- list(X = matrix(1, 2, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb")
    expect_equal(st$log_s, 0)
  })
})

# ==============================================================================
# P.4: .expand_bounds() and .dd_tmb_run_optimizer()
# ==============================================================================

describe(".expand_bounds()", {
  it("returns the default for every name when bounds is NULL", {
    out <- .expand_bounds(NULL, c("beta_k", "beta_k", "log_aux"), -Inf)
    expect_equal(out, c(-Inf, -Inf, -Inf))
  })

  it("applies a named bound to all matching positions", {
    out <- .expand_bounds(c(beta_k = -2), c("beta_k", "beta_k", "log_aux"), -Inf)
    expect_equal(unname(out), c(-2, -2, -Inf))
  })

  it("warns on an unknown bound name", {
    expect_warning(
      .expand_bounds(c(nope = 1), c("beta_k"), Inf),
      "unknown parameter"
    )
  })
})

describe(".dd_tmb_build_map()", {
  it("fixes log_s when !has_s and is NULL when has_s", {
    m1 <- .dd_tmb_build_map(FALSE)
    expect_true(is.list(m1))
    expect_true("log_s" %in% names(m1))
    expect_true(is.factor(m1$log_s))
    expect_true(all(is.na(m1$log_s)))
    expect_null(.dd_tmb_build_map(TRUE))
  })
})

describe(".dd_apply_s_bounds()", {
  it("adds wide log_s bounds only for has_s and only when unset by the user", {
    ctrl <- list(lower = NULL, upper = NULL)
    out <- .dd_apply_s_bounds(ctrl, has_s = TRUE)
    expect_equal(unname(out$lower["log_s"]), log(0.05), tolerance = 1e-12)
    expect_equal(unname(out$upper["log_s"]), log(20), tolerance = 1e-12)
  })

  it("is a no-op when !has_s", {
    ctrl <- list(lower = NULL, upper = NULL)
    out <- .dd_apply_s_bounds(ctrl, has_s = FALSE)
    expect_false("log_s" %in% names(out$lower))
    expect_false("log_s" %in% names(out$upper))
  })

  it("honors a user-supplied log_s bound (user wins)", {
    ctrl <- list(lower = c(log_s = log(0.5)), upper = c(log_s = log(3)))
    out <- .dd_apply_s_bounds(ctrl, has_s = TRUE)
    expect_equal(unname(out$lower["log_s"]), log(0.5), tolerance = 1e-12)
    expect_equal(unname(out$upper["log_s"]), log(3), tolerance = 1e-12)
  })
})

describe(".dd_tmb_run_optimizer()", {
  it("minimizes a quadratic via nlminb and normalizes fields", {
    obj <- list(
      par = c(a = 5, b = -5),
      fn  = function(p) sum((p - c(1, 2))^2),
      gr  = function(p) 2 * (p - c(1, 2))
    )
    ctrl <- list(optimizer = "nlminb", iter_max = 100, eval_max = 200,
                 rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res <- .dd_tmb_run_optimizer(obj, obj$par, ctrl,
                                  user_specified = character(0), verbose = 0)
    expect_equal(unname(res$opt$par), c(1, 2), tolerance = 1e-5)
    expect_equal(res$opt$convergence, 0L)
    expect_true(is.character(res$opt$message))
  })

  it("returns convergence code 99 and Inf objective on optimizer error", {
    obj <- list(par = c(a = 0), fn = function(p) stop("boom"),
                gr  = function(p) p)
    ctrl <- list(optimizer = "nlminb", iter_max = 10, eval_max = 20,
                 rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res <- .dd_tmb_run_optimizer(obj, obj$par, ctrl,
                                  user_specified = character(0), verbose = 0)
    expect_equal(res$opt$convergence, 99L)
    expect_equal(res$opt$objective,   Inf)
  })
})

# ==============================================================================
# P.5: .dd_tmb_multi_start() — 3 starts + phi floor guard
# ==============================================================================

describe(".dd_tmb_multi_start()", {
  it("returns a finite-nll fit with a sane intercept on clean data", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 30, log_k_pop = log(0.02),
                                 sigma_u = 0.5, phi = 10, family = "sltb",
                                 equation = "mazur", seed = 11)
    prep     <- .dd_tmb_prepare_data(sim, "y", "x", "id")
    design   <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "sltb")
    starts   <- .dd_tmb_default_starts(prep, design, "sltb", "mazur")
    ctrl     <- list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
                     rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res      <- .dd_tmb_multi_start(tmb_data, starts, ctrl,
                                    user_specified = character(0), verbose = 0)
    expect_true(is.finite(res$opt$objective))
    beta0   <- res$opt$par[names(res$opt$par) == "beta_k"][1]
    expect_lt(abs(beta0), 20)           # not the k->inf collapse
    log_aux <- res$opt$par[["log_aux"]]
    expect_gt(exp(log_aux), 0.1)        # phi above the floor
  })

  it("keeps log_aux above the phi floor on a boundary-heavy subject", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # A Jarvis-70-style boundary subject (all 0s and 1s) would make the SLT
    # likelihood prefer phi->0, k->inf. The log_aux lower bound (log(0.1))
    # makes that optimum unreachable, so the fit stays sane.
    set.seed(70)
    sim <- simulate_dd_ip(n_subjects = 25, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 8, family = "sltb",
                                 equation = "mazur", seed = 70)
    bad <- data.frame(
      id = factor("boundary"),
      x  = c(7, 30, 180, 365, 730, 1460, 2920),
      y  = c(1,  1,   1,   0,   0,    0,    0)
    )
    sim2 <- rbind(
      data.frame(id = as.character(sim$id), x = sim$x, y = sim$y),
      data.frame(id = as.character(bad$id), x = bad$x, y = bad$y)
    )
    prep     <- .dd_tmb_prepare_data(sim2, "y", "x", "id")
    design   <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "sltb")
    starts   <- .dd_tmb_default_starts(prep, design, "sltb", "mazur")
    ctrl     <- list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
                     rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res      <- .dd_tmb_multi_start(tmb_data, starts, ctrl,
                                    user_specified = character(0), verbose = 0)
    beta0   <- res$opt$par[names(res$opt$par) == "beta_k"][1]
    log_aux <- res$opt$par[["log_aux"]]
    # population k recovered near truth, NOT 414000-style blowup
    expect_lt(exp(beta0), 1)
    # log_aux respects the optimizer lower bound (phi never below the floor)
    expect_gte(exp(log_aux), 0.1 - 1e-8)
  })

  it("RETAINS a genuine low-precision fit (small true phi above the floor)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # True phi = 2 is low precision but well above the 0.1 floor: the bound
    # must NOT discard it. Recovered phi should land near 2, not be pinned to
    # the floor and not be rejected.
    sim <- simulate_dd_ip(n_subjects = 60, log_k_pop = log(0.02),
                                 sigma_u = 0.5, phi = 2, family = "sltb",
                                 equation = "mazur", seed = 202)
    prep     <- .dd_tmb_prepare_data(sim, "y", "x", "id")
    design   <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "sltb")
    starts   <- .dd_tmb_default_starts(prep, design, "sltb", "mazur")
    ctrl     <- list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
                     rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res      <- .dd_tmb_multi_start(tmb_data, starts, ctrl,
                                    user_specified = character(0), verbose = 0)
    phi_hat <- exp(res$opt$par[["log_aux"]])
    # genuine low-phi fit is retained: finite nll, phi in a sane low range and
    # NOT pinned hard to the floor (a recovered ~2 proves it was not discarded)
    expect_true(is.finite(res$opt$objective))
    expect_gt(phi_hat, 0.5)
    expect_lt(phi_hat, 6)
  })
})

# ==============================================================================
# P.6: .dd_tmb_extract_estimates()
# ==============================================================================

describe(".dd_tmb_extract_estimates()", {
  it("returns sane coefficients, renames log_aux, and gates pdHess", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(21)
    sim <- simulate_dd_ip(n_subjects = 30, log_k_pop = log(0.02),
                                 sigma_u = 0.5, phi = 10, family = "sltb",
                                 equation = "mazur", seed = 21)
    prep <- .dd_tmb_prepare_data(sim, "y", "x", "id")
    design <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "sltb")
    starts <- .dd_tmb_default_starts(prep, design, "sltb", "mazur")
    obj <- TMB::MakeADFun(tmb_data, starts, random = "u",
                          map = list(log_s = factor(NA)),
                          DLL = "beezdiscounting", silent = TRUE)
    opt_res <- .dd_tmb_run_optimizer(
      obj, obj$par,
      list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
           rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0),
      character(0), 0)
    est <- .dd_tmb_extract_estimates(obj, opt_res$opt,
                                     n_subjects = prep$n_subjects,
                                     family = "sltb", verbose = 0)
    expect_true("log_phi" %in% names(est$coefficients))
    expect_false("log_aux" %in% names(est$coefficients))
    expect_true("log_phi" %in% names(est$se))
    expect_equal(dim(est$u_hat), c(prep$n_subjects, 1L))
    expect_true(is.logical(est$hessian_pd))
  })

  it("renames log_aux to log_sigma_e under gaussian", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(22)
    sim <- simulate_dd_ip(n_subjects = 25, log_k_pop = log(0.02),
                                 sigma_u = 0.5, sigma_e = 0.08,
                                 family = "gaussian", equation = "mazur",
                                 seed = 22)
    prep <- .dd_tmb_prepare_data(sim, "y", "x", "id")
    design <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "gaussian")
    starts <- .dd_tmb_default_starts(prep, design, "gaussian", "mazur")
    obj <- TMB::MakeADFun(tmb_data, starts, random = "u",
                          map = list(log_s = factor(NA)),
                          DLL = "beezdiscounting", silent = TRUE)
    opt_res <- .dd_tmb_run_optimizer(
      obj, obj$par,
      list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
           rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0),
      character(0), 0)
    est <- .dd_tmb_extract_estimates(obj, opt_res$opt,
                                     n_subjects = prep$n_subjects,
                                     family = "gaussian", verbose = 0)
    expect_true("log_sigma_e" %in% names(est$coefficients))
  })
})

# ==============================================================================
# P.7: .dd_tmb_compute_subject_pars()
# ==============================================================================

describe(".dd_tmb_compute_subject_pars()", {
  it("computes k_i = exp(beta0 + sigma_u * u_i) and omits phi (sltb)", {
    coefs <- c(beta_k = log(0.02), log_sigma_u = log(0.5), log_phi = log(8))
    u_hat <- matrix(c(-1, 0, 2), ncol = 1L)
    # Intercept-only design: one column of 1s, 1 obs per subject.
    design_X <- matrix(1, nrow = 3L, ncol = 1L)
    sp <- .dd_tmb_compute_subject_pars(
      coefficients = coefs, u_hat = u_hat,
      subject_levels = c("a", "b", "c"),
      design_X = design_X, subject_id = c(0L, 1L, 2L),
      equation = "mazur", family = "sltb")
    expect_named(sp, c("id", "u_i", "k"))
    # phi is population-level (MVP): never a subject-level column, even for sltb
    expect_false("phi" %in% names(sp))
    sigma_u <- 0.5
    expect_equal(sp$k, exp(log(0.02) + sigma_u * c(-1, 0, 2)), tolerance = 1e-10)
    expect_equal(sp$u_i, c(-1, 0, 2))
  })

  it("returns id/u_i/k (no phi) for gaussian fits", {
    coefs <- c(beta_k = log(0.02), log_sigma_u = log(0.5), log_sigma_e = log(0.1))
    u_hat <- matrix(c(0, 1), ncol = 1L)
    design_X <- matrix(1, nrow = 2L, ncol = 1L)
    sp <- .dd_tmb_compute_subject_pars(
      coefficients = coefs, u_hat = u_hat,
      subject_levels = c("a", "b"),
      design_X = design_X, subject_id = c(0L, 1L),
      equation = "mazur", family = "gaussian")
    expect_named(sp, c("id", "u_i", "k"))
    expect_false("phi" %in% names(sp))
  })

  it("applies each subject's design ROW (factor-correct k), not the intercept", {
    # Two beta_k columns: intercept + conditionC2 contrast. Subject a is in the
    # reference group (row = c(1, 0)); subject b is in C2 (row = c(1, 1)).
    coefs <- c(beta_k = log(0.02), beta_k = 0.7,
               log_sigma_u = log(0.5), log_phi = log(8))
    u_hat <- matrix(c(0.3, -0.4), ncol = 1L)
    # 2 obs per subject; first row per subject defines the between-subject design.
    design_X <- rbind(c(1, 0), c(1, 0), c(1, 1), c(1, 1))
    subject_id <- c(0L, 0L, 1L, 1L)
    sp <- .dd_tmb_compute_subject_pars(
      coefficients = coefs, u_hat = u_hat,
      subject_levels = c("a", "b"),
      design_X = design_X, subject_id = subject_id,
      equation = "mazur", family = "sltb")
    sigma_u <- 0.5
    # Reference subject: intercept only.
    expect_equal(sp$k[1], exp(log(0.02) + sigma_u * 0.3), tolerance = 1e-12)
    # C2 subject: intercept + conditionC2 contrast (NOT the intercept alone).
    expect_equal(sp$k[2], exp(log(0.02) + 0.7 + sigma_u * -0.4),
                 tolerance = 1e-12)
    # The fix matters: C2 k differs from the intercept-only (buggy) value.
    expect_false(isTRUE(all.equal(
      sp$k[2], exp(log(0.02) + sigma_u * -0.4))))
  })
})

# ==============================================================================
# P.8: fit_dd_tmb() — public API: recovery + object shape + degenerate guard
# ==============================================================================

describe("fit_dd_tmb() recovery", {
  it("recovers population k within 0.15 (sltb x mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(101)
    sim <- simulate_dd_ip(n_subjects = 60, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 10, family = "sltb",
                                 equation = "mazur", seed = 101)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb", verbose = 0)
    expect_s3_class(fit, "beezdiscounting_tmb")
    beta0 <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1])
    truth <- 0.01
    # Explicit relative-error check: |recovered - truth| / truth < TOL.
    # sltb is the principled model; tight 0.15 tolerance is achievable.
    expect_lt(abs(exp(beta0) - truth) / truth, 0.15)
    expect_true(fit$converged)
    expect_true(fit$se_available)
  })

  it("recovers population k within 0.30 (gaussian x mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(102)
    sim <- simulate_dd_ip(n_subjects = 60, log_k_pop = log(0.01),
                                 sigma_u = 0.6, sigma_e = 0.08,
                                 family = "gaussian", equation = "mazur",
                                 seed = 102)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "gaussian", verbose = 0)
    beta0 <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1])
    truth <- 0.01
    # Explicit relative-error check. Gaussian baseline has known [0,1]-clamping
    # bias (observed ~0.19 in practice), so tolerance is wider than sltb.
    expect_lt(abs(exp(beta0) - truth) / truth, 0.30)
    expect_true(fit$converged)
  })

  it("recovers k for the exponential equation (sltb, looser tol)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(103)
    sim <- simulate_dd_ip(n_subjects = 60, log_k_pop = log(0.005),
                                 sigma_u = 0.5, phi = 12, family = "sltb",
                                 equation = "exponential", seed = 103)
    fit <- fit_dd_tmb(sim, equation = "exponential", family = "sltb",
                      verbose = 0)
    expect_true(fit$converged)
    truth_k <- 0.005
    beta0 <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1])
    # Explicit relative-error check. sltb is the principled model; exponential
    # is inherently harder to identify at small k, so tolerance is 0.30.
    expect_lt(abs(exp(beta0) - truth_k) / truth_k, 0.30)
  })

  it("recovers k for gaussian x exponential", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(104)
    sim <- simulate_dd_ip(n_subjects = 60, log_k_pop = log(0.005),
                                 sigma_u = 0.5, sigma_e = 0.06,
                                 family = "gaussian", equation = "exponential",
                                 seed = 104)
    fit <- fit_dd_tmb(sim, equation = "exponential", family = "gaussian",
                      verbose = 0)
    expect_true(fit$converged)
    beta0 <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1])
    truth <- 0.005
    # Explicit relative-error check. Gaussian baseline has known [0,1]-clamping
    # bias (a documented property, not a bug), so tolerance is wider at 0.30.
    expect_lt(abs(exp(beta0) - truth) / truth, 0.30)
  })
})

describe("fit_dd_tmb() object shape", {
  it("assembles the contract fields and subject_pars", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(105)
    sim <- simulate_dd_ip(n_subjects = 40, family = "sltb",
                                 equation = "mazur", seed = 105)
    fit <- fit_dd_tmb(sim, verbose = 0)
    expect_s3_class(fit, "beezdiscounting_tmb")
    expect_true(all(c("call", "opt", "model", "sdr", "param_info",
                      "formula_details", "subject_pars", "loglik", "AIC", "BIC",
                      "converged", "se_available", "data", "data_all",
                      "coercion_info") %in% names(fit)))
    # subject_pars carries NO phi column (phi is population-level in the MVP)
    expect_named(fit$subject_pars, c("id", "u_i", "k"))
    expect_false("phi" %in% names(fit$subject_pars))
    expect_equal(nrow(fit$subject_pars), 40L)
    expect_equal(fit$param_info$family, "sltb")
    expect_equal(fit$param_info$equation, "mazur")
    expect_equal(fit$param_info$n_random_effects, 1L)
    expect_equal(fit$param_info$id_var, "id")
    expect_true("log_phi" %in% names(fit$model$coefficients))
    # call is stored (used by summary() in M.8)
    expect_false(is.null(fit$call))
    expect_true(is.call(fit$call))
    # formula_details stores rhs + contrasts for the emmeans ref grid (E.1)
    expect_true(all(c("X", "rhs", "contrasts") %in% names(fit$formula_details)))
  })

  it("persists factor metadata in param_info and retains design cols in data", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(106)
    sim <- simulate_dd_ip(n_subjects = 40, family = "sltb",
                                 equation = "mazur", seed = 106)
    sim$grp <- factor(rep(c("ctrl", "trt"),
                          length.out = length(unique(sim$id)))[
                            match(sim$id, unique(sim$id))])
    fit <- fit_dd_tmb(sim, factors = "grp", verbose = 0)
    expect_equal(fit$param_info$factors, "grp")
    expect_false(fit$param_info$factor_interaction)
    expect_null(fit$param_info$continuous_covariates)
    # the design factor column survives onto fit$data (row-coherent frame)
    expect_true("grp" %in% names(fit$data))
    expect_equal(nrow(fit$data), nrow(fit$formula_details$X))
  })

  it("keeps a sane population k on a boundary-heavy dataset (regression)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(70)
    sim <- simulate_dd_ip(n_subjects = 25, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 8, family = "sltb",
                                 equation = "mazur", seed = 70)
    bad <- data.frame(
      id = "boundary",
      x = c(7, 30, 180, 365, 730, 1460, 2920),
      y = c(1, 1, 1, 0, 0, 0, 0)
    )
    sim2 <- rbind(
      data.frame(id = as.character(sim$id), x = sim$x, y = sim$y), bad)
    fit <- fit_dd_tmb(sim2, equation = "mazur", family = "sltb", verbose = 0)
    beta0 <- fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1]
    # NOT the k -> 414000 collapse; population k stays plausible
    expect_lt(exp(beta0), 1)
    expect_gt(exp(fit$model$coefficients[["log_phi"]]), 0.1)
  })
})

# ==============================================================================
# R1: se_available requires a positive-definite Hessian
# ==============================================================================

describe("fit_dd_tmb() se_available requires a PD Hessian (R1)", {
  it("a good fit keeps se_available == TRUE with a PD Hessian", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(201)
    sim <- simulate_dd_ip(n_subjects = 50, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 12, family = "sltb",
                                 equation = "mazur", seed = 201)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb", verbose = 0)
    expect_true(isTRUE(fit$hessian_pd))
    expect_true(fit$se_available)
  })

  it("se_available IMPLIES a PD Hessian (never TRUE on a non-PD fit)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(202)
    sim <- simulate_dd_ip(n_subjects = 40, family = "sltb",
                                 equation = "mazur", seed = 202)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb", verbose = 0)
    # The R1 contract: se_available => hessian_pd. (A non-PD Hessian must flip
    # se_available to FALSE so SE-consuming methods don't trust bad SEs.)
    if (isTRUE(fit$se_available)) expect_true(isTRUE(fit$hessian_pd))
    # And it is never TRUE when the sdreport is missing.
    if (is.null(fit$sdr)) expect_false(fit$se_available)
  })
})

# ==============================================================================
# B1: factor/covariate-correct per-subject k (subject_pars uses the design ROW)
# ==============================================================================

describe("fit_dd_tmb() factor-correct subject k (B1)", {
  it("recomputes subject k from the design row + beta + sigma_u*u_i (factor)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # Two conditions with a real shift on log k (C2 discounts faster).
    sim <- simulate_dd_ip(n_subjects = 40, n_conditions = 2,
                                 delta_k = c(0, 0.7), log_k_pop = log(0.01),
                                 sigma_u = 0.5, phi = 12, family = "sltb",
                                 equation = "mazur", seed = 4242)
    fit <- fit_dd_tmb(sim, factors = "condition", verbose = 0)

    coefs   <- fit$model$coefficients
    beta_k  <- unname(coefs[names(coefs) == "beta_k"])
    sigma_u <- exp(unname(coefs[["log_sigma_u"]]))
    Xcn     <- colnames(fit$formula_details$X)
    expect_true(any(grepl("conditionC2", Xcn)))  # design has the C2 contrast
    c2_col  <- which(grepl("conditionC2", Xcn))

    # Map id -> condition from the fitted (row-coherent) frame.
    cond_by_id <- tapply(as.character(fit$data$condition),
                         as.character(fit$data$id),
                         function(z) z[1])
    sp <- fit$subject_pars

    # NON-CIRCULAR: rebuild k_i directly from fit$model$coefficients +
    # fit$formula_details$X (NOT from subject_pars$k), then compare.
    for (i in seq_len(nrow(sp))) {
      this_id   <- as.character(sp$id[i])
      this_cond <- cond_by_id[[this_id]]
      xrow      <- c(1, if (identical(this_cond, "C2")) 1 else 0)
      # Order xrow to the design columns (intercept first, C2 contrast).
      xvec        <- numeric(length(beta_k))
      xvec[1]     <- 1
      xvec[c2_col] <- if (identical(this_cond, "C2")) 1 else 0
      eta_expected <- sum(xvec * beta_k) + sigma_u * sp$u_i[i]
      expect_equal(sp$k[i], exp(eta_expected), tolerance = 1e-6)
    }

    # C1 vs C2 k differ in the direction of beta_conditionC2 (positive => C2
    # has the larger fixed-effect contribution to log k, holding u fixed).
    beta_c2 <- beta_k[c2_col]
    # Compare the fixed-effect (RE = 0) cell k between groups.
    k_c1_fixed <- exp(beta_k[1])
    k_c2_fixed <- exp(beta_k[1] + beta_c2)
    if (beta_c2 > 0) {
      expect_gt(k_c2_fixed, k_c1_fixed)
    } else {
      expect_lt(k_c2_fixed, k_c1_fixed)
    }
  })

  it("predict(type='parameters') and ranef() reflect the corrected k", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 30, n_conditions = 2,
                                 delta_k = c(0, 0.7), log_k_pop = log(0.01),
                                 sigma_u = 0.5, phi = 12, family = "sltb",
                                 equation = "mazur", seed = 99)
    fit <- fit_dd_tmb(sim, factors = "condition", verbose = 0)
    pp  <- predict(fit, type = "parameters")
    rr  <- ranef(fit)
    expect_equal(pp$k, fit$subject_pars$k, tolerance = 1e-12)
    expect_equal(rr$k, fit$subject_pars$k, tolerance = 1e-12)

    # Cross-check against predict(type='response', level='subject'): the per-row
    # subject-conditional k must match subject_pars$k for each subject's rows.
    k_row <- beezdiscounting:::.dd_tmb_predict_k(fit, fit$data, level = "subject")
    k_by_id <- tapply(k_row, as.character(fit$data$id), function(z) z[1])
    sp_k_by_id <- stats::setNames(fit$subject_pars$k,
                                  as.character(fit$subject_pars$id))
    common <- intersect(names(k_by_id), names(sp_k_by_id))
    expect_equal(as.numeric(k_by_id[common]), as.numeric(sp_k_by_id[common]),
                 tolerance = 1e-6)
  })

  it("continuous covariate: subject k reflects the covariate contribution", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 40, log_k_pop = log(0.01),
                                 sigma_u = 0.5, phi = 12, family = "sltb",
                                 equation = "mazur", seed = 77)
    # Between-subject covariate (constant within subject).
    ids <- unique(as.character(sim$id))
    cov_by_id <- stats::setNames(seq_along(ids) / length(ids), ids)
    sim$z <- cov_by_id[as.character(sim$id)]
    fit <- fit_dd_tmb(sim, continuous_covariates = "z", verbose = 0)

    coefs   <- fit$model$coefficients
    beta_k  <- unname(coefs[names(coefs) == "beta_k"])
    sigma_u <- exp(unname(coefs[["log_sigma_u"]]))
    Xcn     <- colnames(fit$formula_details$X)
    z_col   <- which(Xcn == "z")
    z_by_id <- tapply(fit$data$z, as.character(fit$data$id), function(v) v[1])
    sp      <- fit$subject_pars
    for (i in seq_len(nrow(sp))) {
      this_id <- as.character(sp$id[i])
      eta_exp <- beta_k[1] + beta_k[z_col] * z_by_id[[this_id]] +
        sigma_u * sp$u_i[i]
      expect_equal(sp$k[i], exp(eta_exp), tolerance = 1e-6)
    }
  })
})

# ==============================================================================
# B3: phi floor holds on the single optimizer path (multi_start = FALSE)
# ==============================================================================

describe("fit_dd_tmb() phi floor on single path (B3)", {
  it("does not collapse phi to ~0 with multi_start = FALSE on boundary data", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(71)
    sim <- simulate_dd_ip(n_subjects = 25, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 8, family = "sltb",
                                 equation = "mazur", seed = 71)
    bad <- data.frame(
      id = "boundary",
      x = c(7, 30, 180, 365, 730, 1460, 2920),
      y = c(1, 1, 1, 0, 0, 0, 0)
    )
    sim2 <- rbind(
      data.frame(id = as.character(sim$id), x = sim$x, y = sim$y), bad)
    fit <- fit_dd_tmb(sim2, equation = "mazur", family = "sltb",
                      multi_start = FALSE, verbose = 0)
    phi <- exp(fit$model$coefficients[["log_phi"]])
    # The floor (.dd_phi_min = 0.1) must hold on the single path too.
    expect_gte(phi, 0.1 - 1e-6)
  })

  it("the user can override the phi floor via tmb_control$lower$log_aux", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(72)
    sim <- simulate_dd_ip(n_subjects = 30, family = "sltb",
                                 equation = "mazur", seed = 72)
    # A user-supplied log_aux lower bound wins over the default floor; just
    # confirm the fit runs and respects the (looser) override.
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb",
                      multi_start = FALSE,
                      tmb_control = list(lower = c(log_aux = log(0.05))),
                      verbose = 0)
    expect_s3_class(fit, "beezdiscounting_tmb")
    expect_gte(exp(fit$model$coefficients[["log_phi"]]), 0.05 - 1e-6)
  })
})

describe("fit_dd_tmb() log_s mapping (1- vs 2-parameter df)", {
  it("a mazur fit has NO free log_s and unchanged df (3 fixed params)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 40, log_k_pop = log(0.01),
                          sigma_u = 0.5, phi = 12, family = "sltb",
                          equation = "mazur", seed = 301)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb", verbose = 0)
    expect_false("log_s" %in% names(fit$opt$par))
    expect_false("log_s" %in% names(fit$model$coefficients))
    # intercept-only sltb mazur: beta_k(1) + log_sigma_u(1) + log_aux(1) = 3
    expect_equal(length(fit$opt$par), 3L)
    expect_equal(attr(logLik(fit), "df"), 3L)
  })

  it("a green-myerson fit estimates a free log_s (df = 4)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 40, log_k_pop = log(0.01),
                          sigma_u = 0.5, phi = 12, s = 0.7,
                          family = "sltb", equation = "green-myerson",
                          seed = 302)
    fit <- fit_dd_tmb(sim, equation = "green-myerson", family = "sltb",
                      verbose = 0)
    expect_true("log_s" %in% names(fit$opt$par))
    expect_true("log_s" %in% names(fit$model$coefficients))
    expect_equal(length(fit$opt$par), 4L)
    expect_equal(attr(logLik(fit), "df"), 4L)
    expect_true(isTRUE(fit$param_info$has_s))
  })

  it("the start_values$s alias maps to log_s", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 30, s = 1.3, family = "sltb",
                          equation = "rachlin", seed = 303)
    # passing s = 2 must not error (it is converted to log_s = log(2))
    expect_no_error(
      fit <- fit_dd_tmb(sim, equation = "rachlin", family = "sltb",
                        start_values = list(s = 2), verbose = 0)
    )
    expect_true("log_s" %in% names(fit$model$coefficients))
  })

  it("rejects start_values with BOTH s and log_s", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 12, s = 1.3, family = "sltb",
                          equation = "rachlin", seed = 304)
    expect_error(
      fit_dd_tmb(sim, equation = "rachlin", family = "sltb",
                 start_values = list(s = 2, log_s = log(3)), verbose = 0),
      "either 's' or 'log_s'"
    )
  })
})

# ==============================================================================
# 3E: Recovery + reduction-to-Mazur (GM/Rachlin x both families)
# ==============================================================================

describe("fit_dd_tmb() recovers (k, s) for the 2-parameter equations", {
  rich_delays <- c(1, 7, 14, 30, 90, 180, 365, 730, 1460, 2920)

  it("recovers k and s for green-myerson (sltb)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 100, delays = rich_delays,
                          log_k_pop = log(0.01), sigma_u = 0.5, phi = 15,
                          s = 0.6, family = "sltb",
                          equation = "green-myerson", seed = 401)
    fit <- fit_dd_tmb(sim, equation = "green-myerson", family = "sltb",
                      verbose = 0)
    expect_true(fit$converged)
    co <- fit$model$coefficients
    k_hat <- exp(unname(co[names(co) == "beta_k"][1]))
    s_hat <- exp(unname(co[["log_s"]]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.20)
    expect_lt(abs(s_hat - 0.6) / 0.6, 0.25)
  })

  it("recovers k and s for rachlin (sltb)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 100, delays = rich_delays,
                          log_k_pop = log(0.01), sigma_u = 0.5, phi = 15,
                          s = 1.4, family = "sltb",
                          equation = "rachlin", seed = 402)
    fit <- fit_dd_tmb(sim, equation = "rachlin", family = "sltb", verbose = 0)
    expect_true(fit$converged)
    co <- fit$model$coefficients
    k_hat <- exp(unname(co[names(co) == "beta_k"][1]))
    s_hat <- exp(unname(co[["log_s"]]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.20)
    expect_lt(abs(s_hat - 1.4) / 1.4, 0.25)
  })

  it("recovers k and s for green-myerson (gaussian, looser tol)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 100, delays = rich_delays,
                          log_k_pop = log(0.01), sigma_u = 0.5, sigma_e = 0.06,
                          s = 0.6, family = "gaussian",
                          equation = "green-myerson", seed = 403)
    fit <- fit_dd_tmb(sim, equation = "green-myerson", family = "gaussian",
                      verbose = 0)
    expect_true(fit$converged)
    co <- fit$model$coefficients
    k_hat <- exp(unname(co[names(co) == "beta_k"][1]))
    s_hat <- exp(unname(co[["log_s"]]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.30)
    expect_lt(abs(s_hat - 0.6) / 0.6, 0.35)
  })

  it("recovers k and s for rachlin (gaussian, looser tol)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # Rachlin x Gaussian is the noisiest recovery cell: Gaussian-on-bounded-
    # proportions has a systematic [0,1]-clamping bias (k_hat ~14% high), so the
    # k tolerance is the loose 0.30 (as in mazur x gaussian). A 15-seed sweep
    # passes 14/15 with s recovered to <=0.05 rel-err on every seed; this fixed
    # seed sits comfortably inside the band (k rel-err ~0.07, s ~0.02).
    sim <- simulate_dd_ip(n_subjects = 100, delays = rich_delays,
                          log_k_pop = log(0.01), sigma_u = 0.5, sigma_e = 0.06,
                          s = 1.4, family = "gaussian",
                          equation = "rachlin", seed = 407)
    fit <- fit_dd_tmb(sim, equation = "rachlin", family = "gaussian",
                      verbose = 0)
    expect_true(fit$converged)
    co <- fit$model$coefficients
    k_hat <- exp(unname(co[names(co) == "beta_k"][1]))
    s_hat <- exp(unname(co[["log_s"]]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.30)
    expect_lt(abs(s_hat - 1.4) / 1.4, 0.35)
  })
})

describe("fit_dd_tmb() 2-parameter eqns reduce to Mazur at s = 1", {
  it("green-myerson on mazur-simulated data estimates s ~ 1 and recovers k", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 100,
                          delays = c(1, 7, 14, 30, 90, 180, 365, 730, 1460, 2920),
                          log_k_pop = log(0.01), sigma_u = 0.5, phi = 15,
                          family = "sltb", equation = "mazur", seed = 411)
    fit <- fit_dd_tmb(sim, equation = "green-myerson", family = "sltb",
                      verbose = 0)
    expect_true("log_s" %in% names(fit$model$coefficients))  # log_s was free
    s_hat <- exp(unname(fit$model$coefficients[["log_s"]]))
    k_hat <- exp(unname(fit$model$coefficients[
      names(fit$model$coefficients) == "beta_k"][1]))
    expect_lt(abs(s_hat - 1) / 1, 0.20)
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.25)
  })

  it("rachlin on mazur-simulated data estimates s ~ 1 and recovers k", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_ip(n_subjects = 100,
                          delays = c(1, 7, 14, 30, 90, 180, 365, 730, 1460, 2920),
                          log_k_pop = log(0.01), sigma_u = 0.5, phi = 15,
                          family = "sltb", equation = "mazur", seed = 412)
    fit <- fit_dd_tmb(sim, equation = "rachlin", family = "sltb", verbose = 0)
    expect_true("log_s" %in% names(fit$model$coefficients))
    s_hat <- exp(unname(fit$model$coefficients[["log_s"]]))
    k_hat <- exp(unname(fit$model$coefficients[
      names(fit$model$coefficients) == "beta_k"][1]))
    expect_lt(abs(s_hat - 1) / 1, 0.20)
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.25)
  })
})
