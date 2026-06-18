# plot() methods + plot_qq() for the modern modeling tiers. Fit-dependent tests
# are skip_on_cran (matching the TMB suite); brms tests additionally require the
# package and the cached posterior fixtures.

# --- memoized small fits (built once per file run) --------------------------
.plot_cache <- new.env(parent = emptyenv())

tmb_fit <- function() {
  if (is.null(.plot_cache$tmb)) {
    .plot_cache$tmb <- fit_dd_tmb(
      simulate_dd_ip(n_subjects = 12, seed = 1),
      equation = "mazur",
      verbose = 0,
      multi_start = FALSE
    )
  }
  .plot_cache$tmb
}

tmb_factor_fit <- function() {
  if (is.null(.plot_cache$tmbf)) {
    sim <- simulate_dd_ip(
      n_subjects = 24,
      n_conditions = 2,
      delta_k = c(0, 0.6),
      seed = 2
    )
    .plot_cache$tmbf <- fit_dd_tmb(
      sim,
      equation = "mazur",
      factors = "condition",
      verbose = 0,
      multi_start = FALSE
    )
  }
  .plot_cache$tmbf
}

choice_fit <- function() {
  if (is.null(.plot_cache$choice)) {
    .plot_cache$choice <- fit_dd_choice(
      simulate_dd_choice(n_subjects = 14, seed = 1),
      equation = "mazur"
    )
  }
  .plot_cache$choice
}

choice_desc_fit <- function(random_slopes = TRUE) {
  key <- if (random_slopes) "cdesc" else "cdesc_pooled"
  if (is.null(.plot_cache[[key]])) {
    .plot_cache[[key]] <- fit_dd_choice(
      simulate_dd_choice(n_subjects = 14, mode = "descriptive", seed = 3),
      mode = "descriptive",
      random_slopes = random_slopes
    )
  }
  .plot_cache[[key]]
}

tmb_2re_fit <- function() {
  if (is.null(.plot_cache$tmb2re)) {
    sim <- simulate_dd_ip(n_subjects = 24, sigma_phi = 0.3, seed = 5)
    .plot_cache$tmb2re <- fit_dd_tmb(
      sim,
      equation = "mazur",
      random_effects = k + phi ~ 1,
      multi_start = FALSE,
      verbose = 0
    )
  }
  .plot_cache$tmb2re
}

# Number of distinct curve groups in the last (population-line) layer.
n_curve_groups <- function(p) {
  b <- ggplot2::ggplot_build(p)
  length(unique(b$data[[length(p$layers)]]$group))
}

# --- internal helpers -------------------------------------------------------
test_that(".dd_x_seq drops non-positive delays on a log axis", {
  xs <- .dd_x_seq(c(0, 1, 10, 100), n_points = 50, x_trans = "log10")
  expect_true(all(xs > 0))
  expect_length(xs, 50)
  # linear keeps the full span
  xl <- .dd_x_seq(c(0, 100), n_points = 10, x_trans = "linear")
  expect_equal(min(xl), 0)
})

test_that(".dd_drop_nonpos_x removes zero-delay rows on a log axis only", {
  df <- data.frame(x = c(0, 1, 2), y = c(.5, .4, .3))
  expect_equal(nrow(suppressMessages(.dd_drop_nonpos_x(df, "x", "log10"))), 2L)
  expect_equal(nrow(.dd_drop_nonpos_x(df, "x", "linear")), 3L)
})

# --- TMB --------------------------------------------------------------------
test_that("plot.beezdiscounting_tmb returns a ggplot for every type", {
  skip_on_cran()
  fit <- tmb_fit()
  for (ty in c("population", "individual", "parameters", "resid")) {
    p <- plot(fit, type = ty)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("plot.beezdiscounting_tmb: ids subsets the per-subject curves", {
  skip_on_cran()
  fit <- tmb_fit()
  ids <- head(as.character(fit$param_info$subject_levels), 3)
  p <- plot(fit, type = "individual", ids = ids)
  ind <- ggplot2::ggplot_build(p)$data[[2]] # geom: observed(1), ind lines(2)
  expect_equal(length(unique(ind$group)), 3L)
  expect_error(plot(fit, type = "individual", ids = "nope"), "Unknown")
})

test_that("plot.beezdiscounting_tmb: a factor fit draws one curve per level", {
  skip_on_cran()
  fitf <- tmb_factor_fit()
  p <- plot(fitf, type = "population")
  expect_equal(n_curve_groups(p), 2L)
  # `at` conditions to a single level -> a single curve
  p1 <- plot(fitf, type = "population", at = list(condition = "C1"))
  expect_equal(n_curve_groups(p1), 1L)
})

test_that("plot_qq.beezdiscounting_tmb builds from the random intercept", {
  skip_on_cran()
  p <- plot_qq(tmb_fit())
  expect_s3_class(p, "ggplot")
  d <- ggplot2::ggplot_build(p)$data[[1]]
  expect_gt(nrow(d), 0)
})

# --- Choice -----------------------------------------------------------------
test_that("plot.beezdiscounting_choice (structural) covers all types", {
  skip_on_cran()
  fit <- choice_fit()
  for (ty in c("population", "individual", "calibration", "parameters")) {
    expect_no_error(ggplot2::ggplot_build(plot(fit, type = ty)))
  }
  expect_s3_class(plot_qq(fit), "ggplot")
})

test_that("descriptive choice: calibration works, curve types error", {
  skip_on_cran()
  fit <- choice_desc_fit()
  expect_s3_class(plot(fit), "ggplot") # default -> calibration
  expect_no_error(ggplot2::ggplot_build(plot(fit, type = "calibration")))
  expect_error(plot(fit, type = "population"), "descriptive")
  expect_error(plot(fit, type = "parameters"), "descriptive")
})

test_that("plot_qq errors for a pooled descriptive fit (no random effects)", {
  skip_on_cran()
  fit <- choice_desc_fit(random_slopes = FALSE)
  expect_error(plot_qq(fit), "no random effects|No subject random")
})

# --- Comparison forest ------------------------------------------------------
test_that("plot.beezdiscounting_comparison flags CI-excludes-null", {
  skip_on_cran()
  cmp <- get_dd_comparisons(tmb_factor_fit())
  p <- plot(cmp)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_s3_class(plot(cmp, type = "difference"), "ggplot")

  # the flag matches the tidy interval excluding the null
  td <- generics::tidy(cmp, exponentiate = TRUE)
  excl <- (td$conf.low > 1 & td$conf.high > 1) |
    (td$conf.low < 1 & td$conf.high < 1)
  pd <- ggplot2::ggplot_build(p)$data[[3]] # points layer
  expect_equal(sum(excl), sum(pd$colour == "#A25F5F"))
})

# --- brms (cached fixtures) -------------------------------------------------
brms_fixture <- function(name) {
  path <- testthat::test_path("fixtures", "brms", name)
  skip_if_not(file.exists(path), paste("missing fixture", name))
  readRDS(path)
}

test_that("plot.beezdiscounting_brms covers all types (fixture)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  fit <- brms_fixture("fit-mazur-beta.rds")
  for (ty in c("population", "individual", "parameters", "resid")) {
    expect_no_error(ggplot2::ggplot_build(plot(fit, type = ty)))
  }
  # the population curve carries a credible band
  pop <- plot(fit, type = "population")
  has_ribbon <- any(vapply(
    pop$layers,
    function(l) inherits(l$geom, "GeomRibbon"),
    logical(1)
  ))
  expect_true(has_ribbon)
})

test_that("plot.beezdiscounting_choice_brms covers all types (fixture)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  fit <- brms_fixture("fit-choice-mazur-group.rds")
  for (ty in c("population", "individual", "calibration", "parameters")) {
    expect_no_error(ggplot2::ggplot_build(plot(fit, type = ty)))
  }
  # the manually-built posterior implied curve should be monotonically
  # decreasing in delay (a discount curve) and carry a credible band
  cd <- .dd_choice_brms_population_curve(
    fit,
    .dd_curve_newdata(fit, at = NULL, n_points = 50, x_trans = "log10")
  )
  by_grp <- split(cd, cd$.group)
  for (g in by_grp) {
    expect_true(all(diff(g$.value) <= 1e-6))
    expect_true(all(g$.lower <= g$.value & g$.value <= g$.upper))
  }
})

# --- edge cases (Codex review) ----------------------------------------------
test_that("plot_qq.beezdiscounting_tmb facets a 2-RE fit and validates which", {
  skip_on_cran()
  fit <- tmb_2re_fit()
  p <- plot_qq(fit)
  panels <- unique(ggplot2::ggplot_build(p)$data[[1]]$PANEL)
  expect_equal(length(panels), 2L) # re_k + re_phi
  expect_error(plot_qq(fit, which = "not_a_term"), "matched no random-effect")
})

test_that("ids= rejects empty and de-duplicates", {
  skip_on_cran()
  fit <- tmb_fit()
  ids <- as.character(fit$param_info$subject_levels)[1:2]
  p_dup <- plot(fit, type = "individual", ids = c(ids, ids))
  ind <- ggplot2::ggplot_build(p_dup)$data[[2]]
  expect_equal(length(unique(ind$group)), 2L) # deduplicated
  expect_error(plot(fit, type = "individual", ids = character(0)), "empty")
})
