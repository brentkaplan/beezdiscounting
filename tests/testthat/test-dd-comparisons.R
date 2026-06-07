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

# ---------------------------------------------------------------------------
# E.3 get_dd_comparisons() and E.4 tidy.beezdiscounting_comparison()
#
# Contrasts of log k between EMM cells. The discount rate is linear in beta on
# the natural-log scale (log k = X beta_k); contrasts are REPORTED on the log10
# scale (estimate = est_log / log(10)) and, optionally, as multiplicative ratios
# (ratio = exp(est_log)). Every check recomputes the contrast estimate
# INDEPENDENTLY (beta %*% (model.matrix row_i - row_j)) and never trusts the
# engine against itself.
# ---------------------------------------------------------------------------

# Helper: fit a k-condition model and return list(fit, beta, Xc) where Xc is the
# per-condition design (one row per level), independently rebuilt.
.dd_cmp_fit <- function(n_conditions, delta_k, seed) {
  dat <- .simulate_dd_ip_mixed(
    n_subjects = 30 * n_conditions, family = "gaussian", equation = "mazur",
    n_conditions = n_conditions, delta_k = delta_k, seed = seed
  )
  fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                    factors = "condition", multi_start = FALSE, verbose = 0)
  beta <- unname(
    fit$model$coefficients[names(fit$model$coefficients) == "beta_k"]
  )
  lv <- levels(fit$data$condition)
  nd <- data.frame(condition = factor(lv, levels = lv))
  Xc <- stats::model.matrix(~ condition, data = nd,
                            contrasts.arg = fit$formula_details$contrasts)
  Xc <- Xc[, colnames(fit$formula_details$X), drop = FALSE]
  list(fit = fit, beta = beta, Xc = Xc, levels = lv)
}

describe("get_dd_comparisons()", {

  it("returns a classed list(k = list(emmeans, contrasts_log10, contrasts_ratio))", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 101)
    res <- get_dd_comparisons(f$fit)
    expect_s3_class(res, "beezdiscounting_comparison")
    expect_named(res, "k")
    expect_named(res$k, c("emmeans", "contrasts_log10", "contrasts_ratio"))
    expect_named(
      res$k$contrasts_log10,
      c("contrast", "estimate", "std.error", "statistic", "df",
        "conf.low", "conf.high", "p.value")
    )
    expect_named(
      res$k$contrasts_ratio,
      c("contrast", "ratio", "conf.low", "conf.high", "p.value")
    )
    # Metadata attributes.
    expect_equal(attr(res, "backend"), "tmb")
    expect_equal(attr(res, "adjustment_method"), "holm")
    expect_equal(attr(res, "contrast_type_used"), "pairwise")
    expect_equal(attr(res, "contrast_by_used"), "NULL")
  })

  it("pairwise gives choose(n, 2) rows; trt.vs.ctrl gives n - 1", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 102)
    pw <- get_dd_comparisons(f$fit, contrast_type = "pairwise")
    tvc <- get_dd_comparisons(f$fit, contrast_type = "trt.vs.ctrl")
    expect_equal(nrow(pw$k$contrasts_log10), choose(3L, 2L))
    expect_equal(nrow(tvc$k$contrasts_log10), 3L - 1L)
    expect_equal(attr(tvc, "contrast_type_used"), "trt.vs.ctrl")
  })

  it("EMM<->contrast invariance: log10 estimate == (k_log[i] - k_log[j]) / log(10)", {
    f <- .dd_cmp_fit(3, c(0, 0.6, 1.2), seed = 103)
    emm <- get_dd_param_emms(f$fit)            # natural-log k_log per cell
    res <- get_dd_comparisons(f$fit, contrast_type = "pairwise")
    cl <- res$k$contrasts_log10

    # The pairwise order is utils::combn(n, 2) over the EMM cell order. Rebuild
    # the same (i, j) index pairs and compare element-by-element.
    cmb <- utils::combn(nrow(emm), 2L)
    expected_log10 <- (emm$k_log[cmb[1L, ]] - emm$k_log[cmb[2L, ]]) / log(10)
    expect_equal(cl$estimate, expected_log10, tolerance = 1e-8)
  })

  it("recomputes est_log NON-CIRCULARLY from beta %*% (row_i - row_j)", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 104)
    res <- get_dd_comparisons(f$fit, contrast_type = "pairwise")
    cl <- res$k$contrasts_log10

    cmb <- utils::combn(nrow(f$Xc), 2L)
    est_log_truth <- vapply(seq_len(ncol(cmb)), function(k) {
      dx <- f$Xc[cmb[1L, k], ] - f$Xc[cmb[2L, k], ]
      sum(dx * f$beta)
    }, numeric(1))
    # log10 scale.
    expect_equal(cl$estimate, est_log_truth / log(10), tolerance = 1e-8)
    # The ratio block is exp(est_log) on the natural scale.
    expect_equal(res$k$contrasts_ratio$ratio, exp(est_log_truth),
                 tolerance = 1e-8)
  })

  it("non-circular SE: se == sqrt(t(dx) V dx) on the beta_k vcov block", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 105)
    skip_if(is.null(f$fit$sdr) || is.null(f$fit$sdr$cov.fixed))
    res <- get_dd_comparisons(f$fit, contrast_type = "pairwise")
    cl <- res$k$contrasts_log10

    V_full <- as.matrix(f$fit$sdr$cov.fixed)
    bk <- which(names(f$fit$opt$par) == "beta_k")
    V <- V_full[bk, bk, drop = FALSE]
    cmb <- utils::combn(nrow(f$Xc), 2L)
    se_truth <- vapply(seq_len(ncol(cmb)), function(k) {
      dx <- f$Xc[cmb[1L, k], ] - f$Xc[cmb[2L, k], ]
      sqrt(as.numeric(t(dx) %*% V %*% dx))
    }, numeric(1))
    # contrasts_log10 SE is se / log(10).
    expect_equal(cl$std.error, se_truth / log(10), tolerance = 1e-8)
    # statistic is the z on the natural-log scale (scale-invariant ratio).
    expect_equal(cl$statistic, (cl$estimate / cl$std.error), tolerance = 1e-10)
    expect_true(all(is.infinite(cl$df)))
  })

  it("ratio CI brackets the ratio and back-transforms the log CI", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 106)
    res <- get_dd_comparisons(f$fit, ci_level = 0.95)
    cr <- res$k$contrasts_ratio
    expect_true(all(cr$conf.low <= cr$ratio & cr$ratio <= cr$conf.high))
    # ratio CI = exp(log10_CI * log(10)) (same Wald quantile, both scales).
    cl <- res$k$contrasts_log10
    expect_equal(cr$conf.low, exp(cl$conf.low * log(10)), tolerance = 1e-8)
    expect_equal(cr$conf.high, exp(cl$conf.high * log(10)), tolerance = 1e-8)
  })

  it("report_ratios = FALSE omits the contrasts_ratio block", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 107)
    res <- get_dd_comparisons(f$fit, report_ratios = FALSE)
    expect_null(res$k$contrasts_ratio)
    expect_named(res$k, c("emmeans", "contrasts_log10"))
  })

  it("adjust rejects emmeans-only methods and accepts base-R methods", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 108)
    expect_error(get_dd_comparisons(f$fit, adjust = "tukey"),
                 "not a valid p-value adjustment method")
    expect_error(get_dd_comparisons(f$fit, adjust = "sidak"),
                 "not a valid p-value adjustment method")
    # Accepts holm / BH / none.
    expect_s3_class(get_dd_comparisons(f$fit, adjust = "holm"),
                    "beezdiscounting_comparison")
    expect_s3_class(get_dd_comparisons(f$fit, adjust = "BH"),
                    "beezdiscounting_comparison")
    none <- get_dd_comparisons(f$fit, adjust = "none")
    # adjust = "none" leaves raw two-sided z p-values untouched.
    cl <- none$k$contrasts_log10
    expect_equal(cl$p.value, 2 * stats::pnorm(-abs(cl$statistic)),
                 tolerance = 1e-12)
  })

  it("compare_specs picks the retained factor set; unknown names abort", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 109)
    res <- get_dd_comparisons(f$fit, compare_specs = ~ condition)
    expect_equal(nrow(res$k$contrasts_log10), choose(3L, 2L))
    expect_equal(attr(res, "compare_specs_used"), "~condition")
    expect_error(get_dd_comparisons(f$fit, compare_specs = ~ nope),
                 "not in the fit")
    expect_error(get_dd_comparisons(f$fit, compare_specs = "condition"),
                 "one-sided formula")
  })

  it("contrast_by: contrasts restricted per by-cell; p adjusted per cell", {
    # Build a 2-factor interaction model so within-condition contrasts differ
    # across by-cells (site). The simulator emits one factor; attach a second
    # balanced one keyed to subject so it is constant within id.
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 120, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = c(0, 0.5, 1.0), seed = 110
    )
    subj <- levels(dat$id)
    site_lab <- stats::setNames(
      factor(paste0("S", rep_len(1:2, length(subj))), levels = c("S1", "S2")),
      subj
    )
    dat$site <- site_lab[as.character(dat$id)]
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = c("condition", "site"),
                      factor_interaction = TRUE,
                      multi_start = FALSE, verbose = 0)

    res <- get_dd_comparisons(
      fit, compare_specs = ~ condition * site,
      contrast_by = "site", adjust = "holm"
    )
    cl <- res$k$contrasts_log10
    # Within each of 2 sites: choose(3, 2) = 3 condition contrasts => 6 rows.
    expect_equal(nrow(cl), 2L * choose(3L, 2L))
    expect_true("site" %in% names(cl))
    expect_setequal(unique(cl$site), c("S1", "S2"))
    expect_equal(attr(res, "contrast_by_used"), "site")

    # Per-cell p-adjustment: each site's holm adjustment uses only that site's
    # 3 raw p-values. Recompute the S1 block's holm-adjusted p independently and
    # confirm it matches (i.e. S2's p-values did not enter S1's adjustment).
    s1 <- cl[cl$site == "S1", , drop = FALSE]
    raw_s1 <- 2 * stats::pnorm(-abs(s1$statistic))
    expect_equal(s1$p.value, stats::p.adjust(raw_s1, method = "holm"),
                 tolerance = 1e-12)
    # And the GLOBAL holm over all 6 would generally differ -> confirms per-cell.
    raw_all <- 2 * stats::pnorm(-abs(cl$statistic))
    global_holm <- stats::p.adjust(raw_all, method = "holm")
    expect_false(isTRUE(all.equal(cl$p.value, global_holm)))
  })
})

describe("tidy.beezdiscounting_comparison()", {

  it("returns a flat tibble with the cross-backend contract columns", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 120)
    res <- get_dd_comparisons(f$fit)
    td <- generics::tidy(res)
    expect_s3_class(td, "tbl_df")
    expect_named(
      td,
      c("param", "contrast", "estimate", "std.error", "statistic", "df",
        "conf.low", "conf.high", "p.value")
    )
    expect_true(all(td$param == "k"))
    expect_equal(nrow(td), choose(3L, 2L))
    # Flat estimates are the log10 contrasts (default).
    expect_equal(td$estimate, res$k$contrasts_log10$estimate, tolerance = 1e-12)
  })

  it("exponentiate = TRUE presents the ratio scale (10^estimate) with NA SE", {
    f <- .dd_cmp_fit(3, c(0, 0.5, 1.0), seed = 121)
    res <- get_dd_comparisons(f$fit)
    td_log <- generics::tidy(res)
    td_exp <- generics::tidy(res, exponentiate = TRUE)
    expect_equal(td_exp$estimate, 10^td_log$estimate, tolerance = 1e-12)
    expect_equal(td_exp$conf.low, 10^td_log$conf.low, tolerance = 1e-12)
    expect_equal(td_exp$conf.high, 10^td_log$conf.high, tolerance = 1e-12)
    expect_true(all(is.na(td_exp$std.error)))
    # 10^log10 == exp(natural-log) == the ratio block.
    expect_equal(td_exp$estimate, res$k$contrasts_ratio$ratio, tolerance = 1e-8)
  })

  it("carries by-columns through the flat frame when contrast_by is active", {
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 120, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = c(0, 0.5, 1.0), seed = 122
    )
    subj <- levels(dat$id)
    dat$site <- factor(paste0("S", rep_len(1:2, length(subj))),
                       levels = c("S1", "S2"))[
      match(as.character(dat$id), subj)
    ]
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = c("condition", "site"),
                      factor_interaction = TRUE,
                      multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit, compare_specs = ~ condition * site,
                              contrast_by = "site")
    td <- generics::tidy(res)
    expect_true("site" %in% names(td))
    expect_setequal(unique(td$site), c("S1", "S2"))
    expect_equal(nrow(td), 2L * choose(3L, 2L))
  })
})
