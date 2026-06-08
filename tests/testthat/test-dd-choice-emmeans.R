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
