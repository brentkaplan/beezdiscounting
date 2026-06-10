describe("get_dd_param_emms / get_dd_comparisons on a structural choice fit", {
  it("runs and returns per-condition k EMMs (factor on log k)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    skip_if_not_installed("emmeans")
    set.seed(501)
    base <- simulate_dd_choice(n_subjects = 60, log_k_pop = log(0.02),
                               sigma_u = 0.5, gamma = 4, beta0 = 0,
                               equation = "mazur", seed = 501)
    # add a between-subject condition that shifts log k for half the subjects
    ids <- unique(base$id)
    cond <- stats::setNames(rep(c("C1", "C2"), length.out = length(ids)), ids)
    base$condition <- factor(cond[base$id])
    fit <- fit_dd_choice(base, mode = "structural", equation = "mazur",
                         factors = "condition", verbose = 0)
    expect_no_error(emm <- get_dd_param_emms(fit, factors_in_emm = "condition"))
    expect_gte(nrow(as.data.frame(emm)), 2L)
    expect_true(all(c("k", "k_log") %in% names(as.data.frame(emm))))
  })

  it("beta_k covariance block is isolated from log_gamma/beta0", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- simulate_dd_choice(n_subjects = 50, gamma = 4, beta0 = 0.3,
                              equation = "mazur", seed = 502)
    fit <- fit_dd_choice(sim, mode = "structural", intercept = TRUE, verbose = 0)
    # beta_k is first; the cov.fixed block selected by name must be square p x p
    p <- sum(names(fit$opt$par) == "beta_k")
    idx <- which(names(fit$opt$par) == "beta_k")
    expect_identical(idx, seq_len(p))                 # contiguous & first
    if (isTRUE(fit$se_available)) {
      blk <- as.matrix(fit$sdr$cov.fixed)[idx, idx, drop = FALSE]
      expect_equal(dim(blk), c(p, p))
    }
  })
})

describe("emmeans with descriptive params present in the template", {
  it("structural choice still routes beta_k through get_dd_param_emms (factor)", {
    skip_on_cran(); skip_if_not_installed("TMB"); skip_if_not_installed("emmeans")
    dat <- .choice_fit_fixture(n_subjects = 40, seed = 51)
    dat$grp <- factor(rep(c("a", "b"), length.out = nrow(dat)))
    fit <- fit_dd_choice(dat, mode = "structural", equation = "mazur",
                         factors = "grp", verbose = 0)
    expect_equal(which(names(fit$opt$par) == "beta_k")[1], 1L)  # beta_k first
    emm <- get_dd_param_emms(fit, factors_in_emm = "grp")
    expect_gte(nrow(as.data.frame(emm)), 2L)
    expect_true(all(c("k", "k_log") %in% names(as.data.frame(emm))))
  })
  it("descriptive fits are rejected by get_dd_param_emms with a clear message", {
    skip_on_cran(); skip_if_not_installed("TMB")
    dat <- .choice_desc_fixture(n_subjects = 40, seed = 52)
    fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
    expect_error(get_dd_param_emms(fit),
                 regexp = "descriptive|not available|VarCorr|ranef")
  })
  it("get_dd_comparisons also rejects descriptive fits", {
    skip_on_cran(); skip_if_not_installed("TMB")
    dat <- .choice_desc_fixture(n_subjects = 40, seed = 53)
    fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
    expect_error(get_dd_comparisons(fit), regexp = "descriptive|not available|VarCorr|ranef")
  })
})
