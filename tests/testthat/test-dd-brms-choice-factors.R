# TICKET-048: factor designs for fit_dd_choice_brms() + choice EMMs/comparisons --
#
# Mirrors the IP fitter's factor support (build_fixed_rhs design on logk,
# .dd_tmb_build_design guards, real formula_details$X) and removes the v1
# EMM/comparison aborts. Fixture-based tests verify every reported number
# against independent posterior draw summaries (acceptance criterion 2).

skip_if_not_installed("brms")
skip_if_not_installed("posterior")

dd_choice_grp_fixture <- function() {
  path <- testthat::test_path("fixtures", "brms", "fit-choice-mazur-group.rds")
  skip_if_not(file.exists(path), "fixture missing: fit-choice-mazur-group")
  readRDS(path)
}

# --- signature + guards (no sampling) -----------------------------------------------

test_that("fit_dd_choice_brms exposes the factor design arguments", {
  expect_true(all(
    c("factors", "factor_interaction", "continuous_covariates") %in%
      names(formals(fit_dd_choice_brms))
  ))
})

test_that("fit_dd_choice_brms rejects unknown factor columns before sampling", {
  set.seed(9)
  d <- data.frame(
    id = rep(1:6, each = 4), delay = rep(c(1, 7, 30, 90), 6),
    ss_amount = 50, ll_amount = 100, choice = rbinom(24, 1, 0.5)
  )
  expect_error(
    fit_dd_choice_brms(d, factors = "nonexistent", verbose = 0),
    "not found"
  )
})

test_that("fit_dd_choice_brms rejects a rank-deficient design before sampling", {
  set.seed(10)
  d <- data.frame(
    id = rep(1:6, each = 4), delay = rep(c(1, 7, 30, 90), 6),
    ss_amount = 50, ll_amount = 100, choice = rbinom(24, 1, 0.5),
    group = factor(rep(c("a", "a", "b", "b", "b", "b"), each = 4)),
    site = factor(rep(c("x", "x", "x", "x", "y", "y"), each = 4))
  )
  # the a:y cell is empty -> aliased interaction column
  expect_error(
    fit_dd_choice_brms(
      d,
      factors = c("group", "site"), factor_interaction = TRUE, verbose = 0
    ),
    "rank-deficient"
  )
})

# --- fixture-based contracts (factor design) ----------------------------------------

test_that("choice coefficients expand over the factor design under TMB names", {
  fit <- dd_choice_grp_fixture()
  co <- coef(fit)
  expect_identical(sum(names(co) == "beta_k"), 2L)
  expect_true("log_gamma" %in% names(co))
  expect_identical(fit$param_info$factors, "group")
  expect_identical(colnames(fit$formula_details$X), c("(Intercept)", "grouptreat"))
})

test_that("choice tidy expands k rows over the factor design", {
  fit <- dd_choice_grp_fixture()
  t <- tidy(fit)
  expect_identical(
    names(t),
    c("term", "estimate", "std.error", "statistic", "p.value",
      "component", "estimate_scale", "term_display")
  )
  expect_true(all(c("k:(Intercept)", "k:grouptreat", "gamma") %in% t$term))
  expect_identical(t$component[startsWith(t$term, "k:")], rep("fixed", 2L))

  # natural-scale rows are transform-then-summarize over the draws
  draws <- posterior::as_draws_matrix(fit$brmsfit)
  bg <- as.numeric(draws[, "b_logk_grouptreat"])
  expect_equal(
    t$estimate[t$term == "k:grouptreat"], stats::median(exp(bg))
  )
  expect_equal(
    t$std.error[t$term == "k:grouptreat"], stats::sd(exp(bg))
  )
})

test_that("choice confint covers the expanded design with parm aliases", {
  fit <- dd_choice_grp_fixture()
  ci <- confint(fit)
  expect_true(all(c("k:(Intercept)", "k:grouptreat", "gamma") %in% ci$term))
  expect_identical(nrow(confint(fit, parm = "beta_k")), 2L)
  expect_identical(nrow(confint(fit, parm = "k:grouptreat")), 1L)
  expect_identical(nrow(confint(fit, parm = "log_gamma")), 1L)
})

test_that("choice subject_pars apply each subject's design row", {
  fit <- dd_choice_grp_fixture()
  sp <- fit$subject_pars
  expect_identical(nrow(sp), length(fit$param_info$subject_levels))

  draws <- posterior::as_draws_matrix(fit$brmsfit)
  id_treat <- unique(as.character(fit$data$id[fit$data$group == "treat"]))[1]
  lin <- as.numeric(draws[, "b_logk_Intercept"]) +
    as.numeric(draws[, "b_logk_grouptreat"]) +
    as.numeric(draws[, paste0("r_id__logk[", id_treat, ",Intercept]")])
  expect_equal(sp$k[sp$id == id_treat], stats::median(exp(lin)))
})

# --- EMMs / comparisons (v1 aborts removed) -----------------------------------------

test_that("choice EMMs are draws-based over the shared reference grid", {
  fit <- dd_choice_grp_fixture()
  em <- get_dd_param_emms(fit)
  expect_identical(
    names(em),
    c("level", "k", "k_log", "std.error", "conf.low", "conf.high")
  )
  expect_identical(nrow(em), 2L)

  draws <- posterior::as_draws_matrix(fit$brmsfit)
  lin <- as.numeric(draws[, "b_logk_Intercept"]) +
    as.numeric(draws[, "b_logk_grouptreat"])
  i <- grep("treat", em$level)
  expect_equal(em$k[i], stats::median(exp(lin)))
  expect_equal(em$k_log[i], stats::median(lin))
  expect_equal(em$std.error[i], stats::sd(lin))

  # marginalizing everything gives the equal-weight grand mean
  em0 <- get_dd_param_emms(fit, factors_in_emm = character(0))
  expect_identical(nrow(em0), 1L)
  lin0 <- as.numeric(draws[, "b_logk_Intercept"]) +
    0.5 * as.numeric(draws[, "b_logk_grouptreat"])
  expect_equal(em0$k_log, stats::median(lin0))
})

test_that("choice comparisons match independent draw summaries with post.prob", {
  fit <- dd_choice_grp_fixture()
  res <- get_dd_comparisons(fit)
  expect_s3_class(res, "beezdiscounting_comparison")
  expect_identical(attr(res, "backend"), "brms")

  cl <- res$k$contrasts_log10
  expect_identical(nrow(cl), 1L)
  expect_true(is.na(cl$p.value))

  # ctrl - treat contrast draws = -b_grouptreat
  d <- -as.numeric(
    posterior::as_draws_matrix(fit$brmsfit)[, "b_logk_grouptreat"]
  )
  expect_equal(cl$estimate, stats::median(d / log(10)))
  expect_equal(
    cl$post.prob, max(mean(d > 0), mean(d < 0)) + 0.5 * mean(d == 0)
  )
  expect_equal(res$k$contrasts_ratio$ratio, stats::median(exp(d)))

  flat <- tidy(res)
  expect_true("post.prob" %in% names(flat))
})

test_that("explicit adjust requests warn on the choice brms path", {
  fit <- dd_choice_grp_fixture()
  expect_warning(get_dd_comparisons(fit, adjust = "holm"), "posterior")
})
