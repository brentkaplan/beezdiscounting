# Regression backfill for the Phase-1 audit's coverage gaps. Every EMM/predict
# expectation is recomputed INDEPENDENTLY from the fitted coefficients and the
# stored design route (rhs + contrasts), never the function against itself.

skip_on_cran()
skip_if_not_installed("TMB")

# log-k linear predictor for an arbitrary factor/covariate combination, via the
# fit's stored RHS + per-factor contrasts (independent of get_dd_param_emms /
# predict internals).
.lp <- function(fit, newrow) {
  fd <- fit$formula_details
  # Re-level any factor columns to the FITTED levels so model.matrix can apply
  # contrasts on a single-row grid (a bare string would be a 1-level factor).
  for (nm in names(newrow)) {
    if (nm %in% names(fit$data) && is.factor(fit$data[[nm]])) {
      newrow[[nm]] <- factor(as.character(newrow[[nm]]),
                             levels = levels(fit$data[[nm]]))
    }
  }
  X  <- stats::model.matrix(fd$rhs, data = newrow, contrasts.arg = fd$contrasts)
  X  <- X[, colnames(fd$X), drop = FALSE]
  beta <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"])
  as.numeric(X %*% beta)
}

# Balanced full-factorial between-subject design over the given factor levels.
.factorial_data <- function(n_per_cell, levels_list, family = "gaussian", seed = 1) {
  grid    <- expand.grid(levels_list, stringsAsFactors = TRUE)
  n_cells <- nrow(grid)
  sim     <- simulate_dd_ip(n_subjects = n_cells * n_per_cell, family = family,
                            seed = seed)
  ids      <- levels(sim$id)
  cell_idx <- rep(seq_len(n_cells), each = n_per_cell)[seq_along(ids)]
  for (f in names(levels_list)) {
    by_id  <- stats::setNames(as.character(grid[[f]][cell_idx]), ids)
    sim[[f]] <- factor(by_id[as.character(sim$id)], levels = levels(grid[[f]]))
  }
  sim
}

describe("coverage: factors + a continuous covariate together (gap 1)", {
  it("builds a combined design and reconstructs subject k non-circularly", {
    sim <- simulate_dd_ip(n_subjects = 30, n_conditions = 2, delta_k = c(0, 0.6),
                          family = "gaussian", seed = 10)
    ids <- levels(sim$id)
    set.seed(7)
    sim$age <- stats::setNames(stats::rnorm(length(ids)), ids)[as.character(sim$id)]
    fit <- fit_dd_tmb(sim, factors = "condition", continuous_covariates = "age",
                      family = "gaussian", multi_start = FALSE, verbose = 0)
    # the combined design carries BOTH the factor dummy AND the covariate slope
    expect_true(all(c("conditionC2", "age") %in% colnames(fit$formula_details$X)))
    # subject k recomputes from X %*% beta + sigma_u * u_i
    sp      <- fit$subject_pars
    sigma_u <- exp(fit$model$coefficients[["log_sigma_u"]])
    sid     <- match(as.character(fit$data$id), fit$param_info$subject_levels)
    beta    <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"])
    Xd      <- fit$formula_details$X
    k_rec   <- vapply(seq_len(nrow(sp)), function(i) {
      fr <- which(sid == i)[1]
      exp(sum(Xd[fr, ] * beta) + sigma_u * sp$u_i[i])
    }, numeric(1))
    expect_equal(unname(sp$k), unname(k_rec), tolerance = 1e-6)
  })
})

describe("coverage: 3 factors end-to-end (gap 2 / the dropped-3rd-factor regression)", {
  it("places all 3 factors in the design and makes the 3rd EMM-accessible", {
    sim <- .factorial_data(3, list(a = c("a1", "a2"), b = c("b1", "b2"),
                                    c = c("c1", "c2")), seed = 12)
    fit <- fit_dd_tmb(sim, factors = c("a", "b", "c"), family = "gaussian",
                      multi_start = FALSE, verbose = 0)
    # all three factors placed (intercept + one dummy each: <var><level> names)
    expect_equal(ncol(fit$formula_details$X), 4L)
    expect_true(all(c("aa2", "bb2", "cc2") %in% colnames(fit$formula_details$X)))
    em <- get_dd_param_emms(fit, factors_in_emm = "c")   # the previously-dropped 3rd
    expect_equal(nrow(em), 2L)
    # non-circular: EMM(c=ci) marginalizing a,b = exp(mean LP over a,b at c=ci)
    lp_c1 <- mean(c(.lp(fit, data.frame(a = "a1", b = "b1", c = "c1")),
                    .lp(fit, data.frame(a = "a2", b = "b1", c = "c1")),
                    .lp(fit, data.frame(a = "a1", b = "b2", c = "c1")),
                    .lp(fit, data.frame(a = "a2", b = "b2", c = "c1"))))
    expect_equal(em$k[1], exp(lp_c1), tolerance = 1e-6)
  })

  it("fits 3 factors with factor_interaction = TRUE end-to-end (full crossing)", {
    sim <- .factorial_data(4, list(a = c("a1", "a2"), b = c("b1", "b2"),
                                   c = c("c1", "c2")), seed = 13)
    fit <- fit_dd_tmb(sim, factors = c("a", "b", "c"), factor_interaction = TRUE,
                      family = "gaussian", multi_start = FALSE, verbose = 0)
    expect_equal(ncol(fit$formula_details$X), 8L)   # 2^3 full interaction design
    expect_equal(nrow(get_dd_param_emms(fit)), 8L)  # all 8 cells EMM-accessible
  })
})

describe("coverage: at= restricting an OMITTED (marginalized) factor (gap 3)", {
  it("averages only over the at=-selected level, not the full marginal", {
    sim <- .factorial_data(8, list(A = c("a1", "a2"), B = c("b1", "b2")), seed = 21)
    fit <- fit_dd_tmb(sim, factors = c("A", "B"), family = "gaussian",
                      multi_start = FALSE, verbose = 0)
    A_lv <- levels(sim$A)
    emm_at   <- get_dd_param_emms(fit, factors_in_emm = "A", at = list(B = "b1"))
    emm_full <- get_dd_param_emms(fit, factors_in_emm = "A")
    for (i in seq_along(A_lv)) {
      # at B=b1: EMM(A=ai) = exp(LP(ai, b1))
      expect_equal(emm_at$k[i],
                   exp(.lp(fit, data.frame(A = A_lv[i], B = "b1"))),
                   tolerance = 1e-6)
      # full: EMM(A=ai) = exp(mean of LP over both B levels)
      expect_equal(emm_full$k[i],
                   exp(mean(c(.lp(fit, data.frame(A = A_lv[i], B = "b1")),
                              .lp(fit, data.frame(A = A_lv[i], B = "b2"))))),
                   tolerance = 1e-6)
    }
    expect_false(isTRUE(all.equal(emm_at$k, emm_full$k)))   # they genuinely differ
  })
})

describe("coverage: predict() value on a NON-reference factor cell (gap 4)", {
  it("uses the non-reference group's beta (population curve)", {
    sim <- simulate_dd_ip(n_subjects = 24, n_conditions = 2, delta_k = c(0, 0.6),
                          family = "gaussian", seed = 22)
    fit <- fit_dd_tmb(sim, factors = "condition", family = "gaussian",
                      multi_start = FALSE, verbose = 0)
    nd  <- data.frame(x = c(7, 365),
                      condition = factor("C2", levels = c("C1", "C2")))  # non-ref
    p   <- predict(fit, newdata = nd, level = "population")
    k_c2 <- exp(.lp(fit, data.frame(condition = factor("C2", levels = c("C1", "C2")))))
    expect_equal(p$predict.fixed, 1 / (1 + k_c2 * nd$x), tolerance = 1e-6)
  })
})

describe("coverage: trt.vs.ctrl combined with contrast_by (gap 5)", {
  it("computes (nlev-1) per-by-cell contrasts referenced to level 1", {
    sim <- .factorial_data(4, list(condition = c("C1", "C2", "C3"),
                                   site = c("S1", "S2")), seed = 23)
    fit <- fit_dd_tmb(sim, factors = c("condition", "site"), family = "gaussian",
                      multi_start = FALSE, verbose = 0)
    cmp <- get_dd_comparisons(fit, compare_specs = ~ condition + site,
                              contrast_type = "trt.vs.ctrl", contrast_by = "site")
    td  <- generics::tidy(cmp)
    expect_equal(nrow(td), (3L - 1L) * 2L)        # 2 trt-vs-ctrl per site
    expect_true("site" %in% names(td))
    expect_setequal(unique(td$site), c("S1", "S2"))
    # every contrast is "<level> - <reference C1>"
    expect_true(all(grepl("C1$", td$contrast)))
    # per-by-cell p-adjust: within each site, holm is applied over THAT site's
    # raw p-values (2 contrasts), NOT globally over all 4 (non-circular recompute).
    for (s in unique(td$site)) {
      sub   <- td[td$site == s, ]
      raw_p <- 2 * stats::pnorm(-abs(sub$statistic))
      expect_equal(sub$p.value, stats::p.adjust(raw_p, "holm"), tolerance = 1e-8)
    }
  })
})

describe("coverage: unbalanced but full-rank cells (gap 6)", {
  it("equal-weight EMM averaging matches the hand recompute despite unequal n", {
    # 3 subjects in (a1,b1), 9 in (a1,b2), 6 in (a2,b1), 3 in (a2,b2): full-rank,
    # very unbalanced. EMMs use EQUAL weights over cells, not cell sizes.
    counts <- c(a1.b1 = 3, a1.b2 = 9, a2.b1 = 6, a2.b2 = 3)
    sim <- simulate_dd_ip(n_subjects = sum(counts), family = "gaussian", seed = 24)
    ids <- levels(sim$id)
    cell <- rep(names(counts), times = counts)[seq_along(ids)]
    A_by <- stats::setNames(sub("\\..*", "", cell), ids)
    B_by <- stats::setNames(sub(".*\\.", "", cell), ids)
    sim$A <- factor(A_by[as.character(sim$id)])
    sim$B <- factor(B_by[as.character(sim$id)])
    fit <- fit_dd_tmb(sim, factors = c("A", "B"), family = "gaussian",
                      multi_start = FALSE, verbose = 0)
    em <- get_dd_param_emms(fit, factors_in_emm = "A")
    # EMM(A=a1) = exp(mean over B of LP(a1, b)) -- equal weight, ignores n
    expect_equal(em$k[1],
                 exp(mean(c(.lp(fit, data.frame(A = "a1", B = "b1")),
                            .lp(fit, data.frame(A = "a1", B = "b2"))))),
                 tolerance = 1e-6)
  })
})

describe("coverage: boundary-heavy SLT-beta fit, jarvis-independent (gap 7 / spec 4.8)", {
  it("converges with phi >= the floor on simulated boundary-heavy data", {
    # low phi -> indifference points pile up at 0/1 (the SLT-beta motivation),
    # the spec-4.8 degenerate-optimum regime. Simulated (no jarvis dependency).
    sim <- simulate_dd_ip(n_subjects = 40, log_k_pop = log(0.02), sigma_u = 0.8,
                          phi = 1.5, family = "sltb", equation = "mazur", seed = 70)
    expect_gt(mean(sim$y < 0.02 | sim$y > 0.98), 0.1)   # genuinely boundary-heavy
    fit <- fit_dd_tmb(sim, family = "sltb", equation = "mazur", verbose = 0)
    expect_true(fit$converged)
    phi <- exp(fit$model$coefficients[["log_phi"]])
    expect_gte(phi, .dd_phi_min)                    # floor held (no phi -> 0 collapse)
    expect_true(is.finite(exp(fit$model$coefficients[["beta_k"]][1])))
  })
})
