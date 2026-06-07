# Tests for the emmeans foundation of beezdiscounting:
#   .dd_build_emm_ref_grid() (E.1) and get_dd_param_emms() (E.2).
#
# The discount rate is linear in beta on the log scale (log k = X beta_k), so
# every EMM is a deterministic averaging matrix applied to a rebuilt design
# basis, then exp() back-transformed. The checks below recompute the marginal
# means INDEPENDENTLY (beta %*% model.matrix) and never trust the EMM output
# against itself.

skip_on_cran()
skip_if_not_installed("TMB")

describe(".dd_build_emm_ref_grid()", {

  it("returns is_intercept_only for a factor-free, covariate-free fit", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 40, family = "gaussian", equation = "mazur", seed = 1
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      multi_start = FALSE, verbose = 0)
    g <- .dd_build_emm_ref_grid(fit, at = NULL, factors_in_emm = NULL)
    expect_true(g$is_intercept_only)
    expect_null(g$ref_X)
    expect_null(g$level_combos)
  })

  it("builds ref_X = A %*% X_full with equal averaging weights over fitted cols", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 2
    )
    # n_conditions = 2 attaches a `condition` factor (PINNED CONTRACT simulator).
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)

    g <- .dd_build_emm_ref_grid(fit, at = NULL, factors_in_emm = NULL)
    expect_false(g$is_intercept_only)
    # Two observed cells (the two conditions), each a single full-grid row.
    expect_equal(nrow(g$ref_X), 2L)
    expect_equal(ncol(g$ref_X), ncol(fit$formula_details$X))
    expect_setequal(colnames(g$ref_X), colnames(fit$formula_details$X))

    # Non-circular check: rebuild X_full ourselves and confirm ref_X selects
    # the per-cell design rows (A = identity over observed cells when nothing
    # is marginalized).
    full <- expand.grid(condition = levels(fit$data$condition),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    full$condition <- factor(full$condition, levels = levels(fit$data$condition))
    X_full <- stats::model.matrix(~ condition, data = full,
                                  contrasts.arg = attr(fit$formula_details$X, "contrasts"))
    X_full <- X_full[, colnames(fit$formula_details$X), drop = FALSE]
    expect_equal(unname(as.matrix(g$ref_X)), unname(X_full), tolerance = 1e-12)
  })

  it("marginalizes an omitted factor with equal weights (A row sums to 1)", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 3
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    # factors_in_emm = character(0) marginalizes everything to one grand-mean row.
    g <- .dd_build_emm_ref_grid(fit, at = NULL, factors_in_emm = character(0))
    expect_equal(nrow(g$ref_X), 1L)
    # The grand-mean row is the equal-weight average of the per-condition rows.
    full <- expand.grid(condition = levels(fit$data$condition),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    full$condition <- factor(full$condition, levels = levels(fit$data$condition))
    X_full <- stats::model.matrix(~ condition, data = full,
                                  contrasts.arg = attr(fit$formula_details$X, "contrasts"))
    X_full <- X_full[, colnames(fit$formula_details$X), drop = FALSE]
    expect_equal(unname(as.numeric(g$ref_X[1, ])),
                 unname(colMeans(X_full)), tolerance = 1e-12)
  })

  it("pins the FITTED contrasts so a changed options('contrasts') cannot mis-multiply", {
    # Fit with default treatment contrasts.
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 21
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)

    # Flip the session default to sum contrasts; the rebuilt grid must IGNORE
    # this and reproduce the (treatment-coded) fitted column set exactly.
    old <- options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(old), add = TRUE)

    g <- .dd_build_emm_ref_grid(fit, at = NULL, factors_in_emm = NULL)
    expect_equal(colnames(g$ref_X), colnames(fit$formula_details$X))

    # Independent rebuild WITH the stored (fitted) contrasts must match ref_X.
    full <- expand.grid(condition = levels(fit$data$condition),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    full$condition <- factor(full$condition, levels = levels(fit$data$condition))
    X_full <- stats::model.matrix(~ condition, data = full,
                                  contrasts.arg = fit$formula_details$contrasts)
    X_full <- X_full[, colnames(fit$formula_details$X), drop = FALSE]
    expect_equal(unname(as.matrix(g$ref_X)), unname(X_full), tolerance = 1e-12)
  })

  it("aborts when the rebuilt design cannot reproduce the fitted column set", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 22
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)

    # Corrupt the stored fitted column names so the setequal() check fails.
    bad <- fit
    bad$formula_details$X <- bad$formula_details$X[, , drop = FALSE]
    colnames(bad$formula_details$X) <- paste0("zzz_", colnames(bad$formula_details$X))
    expect_error(
      .dd_build_emm_ref_grid(bad, at = NULL, factors_in_emm = NULL),
      "Could not reproduce the fitted design matrix"
    )
  })

  it("validates `at` at the public boundary (.dd_validate_at)", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 23
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    expect_error(
      .dd_build_emm_ref_grid(fit, at = list(condition = "DOES_NOT_EXIST")),
      "not an observed level"
    )
    expect_error(
      .dd_build_emm_ref_grid(fit, at = list(nope = "C1")),
      "Unknown name"
    )
    expect_error(
      .dd_build_emm_ref_grid(fit, at = list("C1")),
      "must be named"
    )
  })
})

describe("get_dd_param_emms()", {

  it("returns the contract columns and back-transforms log k to k = exp(k_log)", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 10
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm <- get_dd_param_emms(fit)
    expect_named(emm, c("level", "k", "k_log", "std.error", "conf.low", "conf.high"))
    expect_equal(nrow(emm), nlevels(fit$data$condition))
    expect_equal(emm$k, exp(emm$k_log), tolerance = 1e-12)
    # CIs back-transform from the log scale (so always positive, bracket k).
    expect_true(all(emm$conf.low > 0))
    expect_true(all(emm$conf.low <= emm$k & emm$k <= emm$conf.high))
  })

  it("recomputes per-cell log-k means NON-CIRCULARLY from coef %*% model.matrix", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 11
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm <- get_dd_param_emms(fit)

    # Truth, by hand: beta_k %*% (per-cell design row), never via emmeans.
    beta <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"])
    lv <- levels(fit$data$condition)
    nd <- data.frame(condition = factor(lv, levels = lv))
    Xc <- stats::model.matrix(~ condition, data = nd,
                              contrasts.arg = attr(fit$formula_details$X, "contrasts"))
    Xc <- Xc[, colnames(fit$formula_details$X), drop = FALSE]
    expected_log_k <- as.numeric(Xc %*% beta)
    expect_equal(emm$k_log, expected_log_k, tolerance = 1e-8)
    # And the back-transform itself, independent of emm$k_log.
    expect_equal(emm$k, exp(expected_log_k), tolerance = 1e-7)
  })

  it("intercept-only fit returns a single (Intercept) row equal to exp(beta_k[1])", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 40, family = "sltb", equation = "mazur", seed = 12
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    emm <- get_dd_param_emms(fit)
    expect_equal(nrow(emm), 1L)
    expect_equal(emm$level, "(Intercept)")
    beta1 <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1])
    expect_equal(emm$k_log, beta1, tolerance = 1e-10)
    expect_equal(emm$k, exp(beta1), tolerance = 1e-10)
  })

  it("marginalizes an omitted factor: factors_in_emm = character(0) gives one cell", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 13
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm0 <- get_dd_param_emms(fit, factors_in_emm = character(0))
    expect_equal(nrow(emm0), 1L)
    # The grand mean of the per-cell log-k EMMs (equal weights).
    emm <- get_dd_param_emms(fit)
    expect_equal(emm0$k_log, mean(emm$k_log), tolerance = 1e-8)
  })

  it("conditions on a continuous covariate via `at` (exp(beta_int + beta_cov*value))", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 80, family = "gaussian", equation = "mazur", seed = 30
    )
    # Attach a deterministic continuous covariate (constant within subject).
    set.seed(31)
    subj_cov <- stats::setNames(stats::rnorm(nlevels(dat$id)), levels(dat$id))
    dat$age <- as.numeric(subj_cov[as.character(dat$id)])

    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      continuous_covariates = "age",
                      multi_start = FALSE, verbose = 0)

    beta <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"])
    cols <- colnames(fit$formula_details$X)
    b_int <- beta[match("(Intercept)", cols)]
    b_age <- beta[match("age", cols)]

    val1 <- 1.0
    val2 <- 2.0
    e1 <- get_dd_param_emms(fit, at = list(age = val1))
    e2 <- get_dd_param_emms(fit, at = list(age = val2))

    expect_equal(nrow(e1), 1L)
    expect_equal(e1$k, exp(b_int + b_age * val1), tolerance = 1e-6)
    expect_equal(e2$k, exp(b_int + b_age * val2), tolerance = 1e-6)
    # Monotone in the covariate value (sign follows b_age).
    if (b_age >= 0) {
      expect_gt(e2$k_log, e1$k_log)
    } else {
      expect_lt(e2$k_log, e1$k_log)
    }
    # Label reflects the covariate=value.
    expect_match(e1$level, "^age=")
  })

  it("uses the beta_k block of sdr\\$cov.fixed for the SE (not the full vcov)", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = c(0, 0.8), seed = 14
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    skip_if(is.null(fit$sdr) || is.null(fit$sdr$cov.fixed))

    emm <- get_dd_param_emms(fit)

    # Independent SE recompute: sub-index cov.fixed to the beta_k block by
    # opt$par naming, then se_x = sqrt(t(x) V x) for each ref row.
    V_full <- as.matrix(fit$sdr$cov.fixed)
    bk <- which(names(fit$opt$par) == "beta_k")
    V <- V_full[bk, bk, drop = FALSE]
    g <- .dd_build_emm_ref_grid(fit, factors_in_emm = NULL, validate = FALSE)
    ref_X <- g$ref_X
    se_expected <- sqrt(diag(ref_X %*% V %*% t(ref_X)))
    expect_equal(emm$std.error, as.numeric(se_expected), tolerance = 1e-10)
  })
})
