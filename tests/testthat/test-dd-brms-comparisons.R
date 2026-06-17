# TICKET-041: draws-based EMMs/comparisons + init = "tmb" for the dd brms tier --

skip_if_not_installed("brms")
skip_on_ci()  # brms fits real Stan models; too slow/fragile under covr on CI (run locally)
skip_if_not_installed("posterior")

dd_grp_fixture <- function() {
  path <- testthat::test_path("fixtures", "brms", "fit-mazur-beta-group.rds")
  skip_if_not(file.exists(path), "fixture missing: fit-mazur-beta-group")
  readRDS(path)
}

test_that("get_dd_param_emms is draws-based for brms fits", {
  fit <- dd_grp_fixture()
  em <- get_dd_param_emms(fit)

  expect_identical(
    names(em),
    c("level", "k", "k_log", "std.error", "conf.low", "conf.high")
  )
  expect_identical(nrow(em), 2L)

  # behavioral: treat-level k = median of exp(intercept + grouptreat) draws
  draws <- posterior::as_draws_matrix(fit$brmsfit)
  lin <- as.numeric(draws[, "b_logk_Intercept"]) +
    as.numeric(draws[, "b_logk_grouptreat"])
  i <- grep("treat", em$level)
  expect_equal(em$k[i], stats::median(exp(lin)))
  expect_equal(em$k_log[i], stats::median(lin))
  expect_equal(em$std.error[i], stats::sd(lin))
})

test_that("get_dd_comparisons matches independent draw summaries with post.prob", {
  fit <- dd_grp_fixture()
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
  expect_equal(cl$post.prob, max(mean(d > 0), mean(d < 0)))
  expect_equal(res$k$contrasts_ratio$ratio, stats::median(exp(d)))

  # tidy carries post.prob through the flat schema
  flat <- tidy(res)
  expect_true("post.prob" %in% names(flat))
})

test_that("explicit adjust requests warn on the brms path", {
  fit <- dd_grp_fixture()
  expect_warning(get_dd_comparisons(fit, adjust = "holm"), "posterior")
})

test_that("intercept-only brms fits return empty contrasts", {
  path <- testthat::test_path("fixtures", "brms", "fit-mazur-beta.rds")
  skip_if_not(file.exists(path), "fixture missing")
  fit <- readRDS(path)
  res <- get_dd_comparisons(fit)
  expect_identical(nrow(res$k$contrasts_log10), 0L)
  em <- get_dd_param_emms(fit)
  expect_identical(em$level, "(Intercept)")
})

test_that("redundant contrast_by is ignored with global contrasts (Codex 041-B1)", {
  fit <- dd_grp_fixture()
  expect_message(
    res <- get_dd_comparisons(fit, contrast_by = "group"),
    "redundant"
  )
  cl <- res$k$contrasts_log10
  expect_identical(nrow(cl), 1L) # global group contrast, not empty
  expect_false("group" %in% names(cl))
  expect_identical(attr(res, "contrast_by_used"), "NULL")
})

test_that("post.prob is tie-aware (Codex 041-R1)", {
  pp <- beezdiscounting:::.dd_brms_post_prob
  expect_identical(pp(rep(0, 8)), 0.5)
  expect_identical(pp(c(1, 1, 1, -1)), 0.75)
  expect_identical(pp(c(1, 1, 0, 0)), 0.75)
})

test_that("init = 'tmb' uses the FULL beta_k vector for factor designs (Codex 041-R2)", {
  fit <- dd_grp_fixture()
  d <- fit$data

  spec <- beezdiscounting:::.dd_brms_formula(
    "mazur", "beta",
    factors = "group", data = d
  )
  inits <- beezdiscounting:::.dd_brms_build_inits(
    init = "tmb", spec = spec, data = d, chains = 1, seed = 1,
    autoscale_info = NULL, family = "beta",
    factors = "group", factor_interaction = FALSE,
    continuous_covariates = NULL, equation = "mazur"
  )
  expect_length(inits[[1]]$b_logk, 2L)

  tmb <- fit_dd_tmb(
    d,
    equation = "mazur", family = "sltb", factors = "group", verbose = 0
  )
  beta_k <- unname(tmb$model$coefficients[
    names(tmb$model$coefficients) == "beta_k"
  ])
  # BOTH coefficients centered at the TMB MLE (jitter sd 0.1), not just
  # the intercept with zeros
  expect_lt(max(abs(as.numeric(inits[[1]]$b_logk) - beta_k)), 0.5)
  expect_gt(abs(beta_k[2]), 0.5) # the group effect is genuinely nonzero
})

test_that("init = 'tmb' centers the inits at a TMB pre-fit (no sampling)", {
  set.seed(31)
  delays <- c(1, 7, 30, 90, 180, 365)
  d <- expand.grid(id = factor(1:8), x = delays)
  k_i <- exp(log(0.02) + rnorm(8, 0, 0.4))
  d$y <- pmin(pmax(1 / (1 + k_i[d$id] * d$x) + rnorm(48, 0, 0.05), 0.01), 0.99)

  spec <- beezdiscounting:::.dd_brms_formula("mazur", "beta")
  inits <- beezdiscounting:::.dd_brms_build_inits(
    init = "tmb", spec = spec, data = d, chains = 2, seed = 1,
    autoscale_info = NULL, family = "beta",
    factors = NULL, factor_interaction = FALSE,
    continuous_covariates = NULL, equation = "mazur"
  )
  expect_length(inits, 2L)
  expect_true(all(c("b_logk", "sd_1", "z_1", "phi") %in% names(inits[[1]])))

  # centers track the TMB sltb MLE (jitter sd = 0.1)
  tmb <- fit_dd_tmb(d, equation = "mazur", family = "sltb", verbose = 0)
  expect_lt(
    abs(as.numeric(inits[[1]]$b_logk[1]) - tmb$model$coefficients[["beta_k"]]),
    0.5
  )
})
