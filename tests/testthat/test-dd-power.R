# Heavy file: runs only with BEEZ_FULL_TESTS=true (full-tests.yaml tri-OS job
# and a local `BEEZ_FULL_TESTS=true` run); skipped in every ungated run
# (test-coverage.yaml, plain devtools::test()) -- under covr's instrumented
# TMB build the thousands of refits below take > 6 h. The no-model-fit tests
# live in test-dd-power-fast.R. See tests/testthat/helper-full-tests.R.
.skip_unless_full_tests()

# =============================================================================
# Monte Carlo power analysis: power_discounting() / find_n_discounting()
#
# Mirrors beezdemand's power_demand() validity battery (between-subject
# design, delta_k on log k). See vignette("power-analysis") for the
# "Validity and Limitations" statement backed by these tests.
#
#   1. Type I error calibration (the load-bearing test). n_sim is computed
#      from the tolerance, not guessed: to assert |rate - .05| <= .02 with a
#      >= 3-sigma acceptance band under correct calibration, n_sim >=
#      9 * .05 * .95 / .02^2 = 1069 -> n_sim = 1200 (band [.03, .07], which
#      excludes both 0.5x and 1.5x the nominal rate). Seeds and configs are
#      fixed BEFORE the first run; tolerances are not adjusted afterwards.
#      The engine refers the Wald statistic to t(n_subjects - 2) by default
#      (two-sample df): the beezdemand calibration battery showed the
#      asymptotic z-test is anticonservative at study-relevant N, and the
#      same correction is mirrored here; this test validates it for this
#      package's between-subject design.
#   2. Convergence handling: unusable replicates are excluded from the
#      denominator and surfaced, never counted as misses.
#   3. Closed-form benchmark vs pwr::pwr.t.test() (even N so the round-robin
#      condition assignment gives equal groups).
#   4. Monotonicity sweeps (one fixed direction per factor).
#   5. Reproducibility: identical seed -> identical result.
# =============================================================================


test_that("the RNG state is left exactly as found", {
  skip_on_cran()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 0.5),
    n_sim = 1,
    seed = 5,
    verbose = FALSE
  ))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  set.seed(999)
  before <- get(".Random.seed", envir = globalenv())
  invisible(power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 0.5),
    n_sim = 1,
    seed = 5,
    verbose = FALSE
  ))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

# -----------------------------------------------------------------------------
# Return structure + p-value/CI verdict agreement
# -----------------------------------------------------------------------------

test_that("power_discounting returns the documented structure", {
  skip_on_cran()
  res <- power_discounting(
    n_subjects = 12,
    effect = list(delta_k = 1),
    n_sim = 8,
    seed = 101,
    verbose = FALSE
  )

  expect_s3_class(res, "beezdiscounting_power")
  expect_true(all(
    c(
      "power",
      "power_mc_ci",
      "hit_rate_p",
      "hit_rate_ci",
      "n_sim",
      "n_converged",
      "n_hessian_pd",
      "n_used",
      "alpha",
      "effect",
      "target_term",
      "design",
      "n_subjects",
      "replicates"
    ) %in%
      names(res)
  ))
  expect_equal(res$n_sim, 8)
  expect_equal(nrow(res$replicates), 8)
  expect_equal(res$target_term, "k:conditionC2")
  expect_gte(res$power, 0)
  expect_lte(res$power, 1)
  expect_lte(res$power_mc_ci[1], res$power)
  expect_gte(res$power_mc_ci[2], res$power)

  used <- res$replicates[res$replicates$status == "ok", ]
  expect_equal(used$hit_p, used$hit_ci)

  expect_output(print(res), "Monte Carlo power")
})

# -----------------------------------------------------------------------------
# Estimate recovery on the simulator's log-k scale
# -----------------------------------------------------------------------------

test_that("power_discounting extracts delta_k on the simulator's log scale", {
  skip_on_cran()
  delta <- 0.8
  res <- power_discounting(
    n_subjects = 40,
    effect = list(delta_k = delta),
    n_sim = 12,
    seed = 103,
    verbose = FALSE
  )
  used <- res$replicates[res$replicates$status == "ok", ]
  expect_gt(nrow(used), 6)
  # If extraction were on the exp scale this would sit near exp(0.8) = 2.23.
  expect_equal(mean(used$estimate), delta, tolerance = 0.2)
})

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------

test_that("power_discounting is exactly reproducible under a seed", {
  skip_on_cran()
  res1 <- power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 0.8),
    n_sim = 5,
    seed = 42,
    verbose = FALSE
  )
  res2 <- power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 0.8),
    n_sim = 5,
    seed = 42,
    verbose = FALSE
  )
  expect_identical(res1$replicates, res2$replicates)
  expect_identical(res1$power, res2$power)

  res3 <- power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 0.8),
    n_sim = 5,
    seed = 43,
    verbose = FALSE
  )
  expect_false(identical(res1$replicates, res3$replicates))
})

# -----------------------------------------------------------------------------
# Convergence handling, tested not assumed
# -----------------------------------------------------------------------------

test_that("all-nonconverged replicates yield NA power and n_converged = 0", {
  skip_on_cran()
  res <- suppressWarnings(power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 0.5),
    n_sim = 4,
    seed = 44,
    verbose = FALSE,
    tmb_control = list(iter_max = 1, eval_max = 2)
  ))
  expect_equal(res$n_converged, 0)
  expect_equal(res$n_used, 0)
  expect_true(is.na(res$power))
  expect_true(all(is.na(res$replicates$hit_p)))
  # Nonconverged replicates surface the optimizer's message and keep the
  # "nonconverged" status (they are not execution errors).
  expect_equal(unique(res$replicates$status), "nonconverged")
  expect_true(all(!is.na(res$replicates$message)))
})

test_that("partially non-converged runs exclude failures from the denominator", {
  skip_on_cran()
  # Config found empirically to give a stable mix of usable and failed fits
  # under this seed (tiny N per group, extreme effect, huge heterogeneity,
  # low-precision SLT-beta, two delays).
  res <- suppressWarnings(power_discounting(
    n_subjects = 4,
    effect = list(delta_k = 6),
    design = list(sigma_u = 2.5, phi = 1.5, delays = c(7, 30)),
    n_sim = 30,
    seed = 20260720,
    verbose = FALSE
  ))
  expect_lt(res$n_used, res$n_sim)
  expect_gt(res$n_used, 0)
  used <- res$replicates[res$replicates$status == "ok", ]
  expect_equal(nrow(used), res$n_used)
  expect_equal(res$power, mean(used$hit_ci))
})

# -----------------------------------------------------------------------------
# Type I error calibration -- THE load-bearing test (preregistered).
# -----------------------------------------------------------------------------

test_that("Type I error is calibrated at nominal alpha under the null", {
  skip_on_cran()
  res <- power_discounting(
    n_subjects = 30,
    effect = list(delta_k = 0),
    n_sim = 1200,
    seed = 20260721,
    verbose = FALSE
  )
  expect_gte(res$n_used / res$n_sim, 0.95)
  expect_gte(res$hit_rate_p, 0.03)
  expect_lte(res$hit_rate_p, 0.07)
  expect_gte(res$hit_rate_ci, 0.03)
  expect_lte(res$hit_rate_ci, 0.07)
})

test_that("Type I error holds at a realistic N", {
  skip_on_cran()
  # Band: .05 +/- 3.18 * sqrt(.05 * .95 / 400) = [.015, .085] (preregistered).
  res <- power_discounting(
    n_subjects = 60,
    effect = list(delta_k = 0),
    n_sim = 400,
    seed = 20260722,
    verbose = FALSE
  )
  expect_gte(res$n_used / res$n_sim, 0.95)
  expect_gte(res$hit_rate_p, 0.015)
  expect_lte(res$hit_rate_p, 0.085)
})

# -----------------------------------------------------------------------------
# Closed-form benchmark. With family = "gaussian", tiny residual error, and
# many delays, each subject's log k is recovered nearly exactly, so the
# between-subject delta_k test reduces to a two-sample t-test on log k with
# d = delta_k / sigma_u and equal groups (even N; round-robin assignment).
# Preregistered: analytic power computed BEFORE the MC run; tolerance 0.10 =
# 3 MC SDs at n_sim = 400 (~.065) plus z-vs-t slack at n = 20/group.
# -----------------------------------------------------------------------------

test_that("Monte Carlo power matches analytic power in a degenerate design", {
  skip_on_cran()
  skip_if_not_installed("pwr")

  n_subj <- 40 # even: round-robin gives 20 per condition
  sigma_u <- 0.5
  delta <- 0.45
  analytic <- pwr::pwr.t.test(
    n = n_subj / 2,
    d = delta / sigma_u,
    sig.level = 0.05,
    type = "two.sample"
  )$power

  res <- power_discounting(
    n_subjects = n_subj,
    effect = list(delta_k = delta),
    design = list(sigma_u = sigma_u, sigma_e = 0.02),
    family = "gaussian",
    n_sim = 400,
    seed = 20260723,
    verbose = FALSE
  )
  expect_gte(res$n_used / res$n_sim, 0.95)
  expect_equal(res$power, analytic, tolerance = 0.10)
})

# -----------------------------------------------------------------------------
# Monotonicity sanity sweeps
# -----------------------------------------------------------------------------

test_that("power increases with n_subjects", {
  skip_on_cran()
  p_small <- power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 0.5),
    n_sim = 100,
    seed = 501,
    verbose = FALSE
  )$power
  p_large <- power_discounting(
    n_subjects = 40,
    effect = list(delta_k = 0.5),
    n_sim = 100,
    seed = 502,
    verbose = FALSE
  )$power
  expect_gt(p_large, p_small)
})

test_that("power increases with effect size", {
  skip_on_cran()
  p_small <- power_discounting(
    n_subjects = 20,
    effect = list(delta_k = 0.3),
    n_sim = 100,
    seed = 503,
    verbose = FALSE
  )$power
  p_large <- power_discounting(
    n_subjects = 20,
    effect = list(delta_k = 1.2),
    n_sim = 100,
    seed = 504,
    verbose = FALSE
  )$power
  expect_gt(p_large, p_small)
})

test_that("power decreases as sigma_u grows", {
  skip_on_cran()
  p_low_var <- power_discounting(
    n_subjects = 20,
    effect = list(delta_k = 0.5),
    design = list(sigma_u = 0.3),
    n_sim = 100,
    seed = 505,
    verbose = FALSE
  )$power
  p_high_var <- power_discounting(
    n_subjects = 20,
    effect = list(delta_k = 0.5),
    design = list(sigma_u = 1.2),
    n_sim = 100,
    seed = 506,
    verbose = FALSE
  )$power
  expect_lt(p_high_var, p_low_var)
})

# -----------------------------------------------------------------------------
# find_n_discounting
# -----------------------------------------------------------------------------

test_that("find_n_discounting finds a plausible minimum N for a large effect", {
  skip_on_cran()
  res <- find_n_discounting(
    target_power = 0.8,
    effect = list(delta_k = 2),
    n_range = c(4, 40),
    n_sim = 60,
    seed = 401,
    verbose = FALSE
  )
  expect_s3_class(res, "beezdiscounting_power_n")
  expect_true(all(
    c(
      "n",
      "target_power",
      "status",
      "uncertain",
      "evaluations"
    ) %in%
      names(res)
  ))
  expect_gte(res$n, 4)
  expect_lte(res$n, 40)
  sel <- res$evaluations[res$evaluations$n_subjects == res$n, ]
  expect_gte(max(sel$power), 0.7)
  # The final (confirmation) look at the selected N must have cleared the
  # target -- decisively, or by point estimate when flagged uncertain.
  expect_true(sel$decision[nrow(sel)] %in% c("above", "ambiguous_above"))
  if (!res$uncertain) expect_equal(sel$decision[nrow(sel)], "above")
  expect_output(print(res), "Monte Carlo uncertainty|minimum")
})

test_that("find_n_discounting reports an unreachable target instead of extrapolating", {
  skip_on_cran()
  expect_error(
    find_n_discounting(
      target_power = 0.95,
      effect = list(delta_k = 0.05),
      n_range = c(4, 6),
      n_sim = 40,
      seed = 402,
      verbose = FALSE
    ),
    "not reach|unreachable"
  )
})

test_that("target extraction is robust to a global contr.sum option", {
  skip_on_cran()
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  res <- power_discounting(
    n_subjects = 10,
    effect = list(delta_k = 1),
    n_sim = 2,
    seed = 7,
    verbose = FALSE
  )
  # The engine forces treatment coding via a scoped options override, so the
  # target term is found and the estimand is still C2 - C1 on the log-k scale.
  expect_equal(res$n_used, 2)
  expect_equal(sum(is.na(res$replicates$estimate)), 0)
})
