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

# -----------------------------------------------------------------------------
# Input validation (cheap, no MC)
# -----------------------------------------------------------------------------

test_that("power_discounting validates the effect specification", {
  expect_error(
    power_discounting(n_subjects = 10, effect = list(delta_k = NULL)),
    "exactly one"
  )
  expect_error(
    power_discounting(n_subjects = 10, effect = list(bad_name = 0.5)),
    "delta_k"
  )
  expect_error(
    power_discounting(n_subjects = 10, effect = list(delta_k = "big")),
    "single finite number"
  )
})

test_that("power_discounting validates design elements", {
  expect_error(
    power_discounting(
      n_subjects = 10,
      effect = list(delta_k = 0.5),
      design = list(not_a_design_arg = 1)
    ),
    "not_a_design_arg"
  )
  expect_error(
    power_discounting(
      n_subjects = 10,
      effect = list(delta_k = 0.5),
      design = list(sigma_u = -1)
    ),
    "sigma_u"
  )
})

test_that("power_discounting restricts equation to the v1 scope", {
  expect_error(
    power_discounting(
      n_subjects = 10,
      effect = list(delta_k = 0.5),
      equation = "rachlin"
    )
  )
})

test_that("power_discounting rejects unnamed or duplicated list elements", {
  expect_error(
    power_discounting(n_subjects = 10, effect = list(0.5)),
    "named"
  )
  expect_error(
    power_discounting(
      n_subjects = 10,
      effect = list(delta_k = 0.5, delta_k = 9)
    ),
    "duplicated"
  )
  expect_error(
    power_discounting(
      n_subjects = 10,
      effect = list(delta_k = 0.5),
      design = list(0.2)
    ),
    "named"
  )
})

test_that("power_discounting rejects designs whose default df is not positive", {
  expect_error(
    power_discounting(n_subjects = 2, effect = list(delta_k = 0.5)),
    "df"
  )
  expect_error(
    find_n_discounting(
      target_power = 0.8,
      effect = list(delta_k = 0.5),
      n_range = c(2, 10)
    ),
    "n_range"
  )
})

test_that("power_discounting validates the seed range", {
  expect_error(
    power_discounting(
      n_subjects = 10,
      effect = list(delta_k = 0.5),
      seed = 2^31
    ),
    "seed"
  )
})

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

test_that("power_discounting validates scalar arguments", {
  expect_error(
    power_discounting(n_subjects = 1, effect = list(delta_k = 0.5)),
    "n_subjects"
  )
  expect_error(
    power_discounting(
      n_subjects = 10,
      effect = list(delta_k = 0.5),
      alpha = 1.2
    ),
    "alpha"
  )
  expect_error(
    power_discounting(n_subjects = 10, effect = list(delta_k = 0.5), n_sim = 0),
    "n_sim"
  )
})

# -----------------------------------------------------------------------------
# Wilson interval helper (pure, unit-level)
# -----------------------------------------------------------------------------

test_that(".dd_power_wilson_ci matches known values and handles edges", {
  ci <- beezdiscounting:::.dd_power_wilson_ci(5, 10)
  expect_equal(ci, c(0.2366, 0.7634), tolerance = 1e-3)

  ci0 <- beezdiscounting:::.dd_power_wilson_ci(0, 10)
  expect_equal(ci0[1], 0)
  expect_equal(ci0[2], 0.2775, tolerance = 1e-3)

  ci_all <- beezdiscounting:::.dd_power_wilson_ci(10, 10)
  expect_equal(ci_all[2], 1)
  expect_gt(ci_all[1], 0.7)

  ci_none <- beezdiscounting:::.dd_power_wilson_ci(0, 0)
  expect_equal(ci_none, c(NA_real_, NA_real_))
})

# -----------------------------------------------------------------------------
# Replicate classification (pure, unit-level)
# -----------------------------------------------------------------------------

test_that(".dd_power_rep_row classifies unusable fits and never counts them as misses", {
  ok <- beezdiscounting:::.dd_power_rep_row(
    converged = TRUE,
    hessian_pd = TRUE,
    estimate = 0.5,
    se = 0.1,
    alpha = 0.05
  )
  expect_equal(ok$status, "ok")
  expect_equal(ok$hit_p, TRUE)
  expect_equal(ok$hit_ci, TRUE)

  nonconv <- beezdiscounting:::.dd_power_rep_row(
    converged = FALSE,
    hessian_pd = FALSE,
    estimate = 0.5,
    se = 0.1,
    alpha = 0.05
  )
  expect_equal(nonconv$status, "nonconverged")
  expect_equal(nonconv$hit_p, NA)

  bad_hess <- beezdiscounting:::.dd_power_rep_row(
    converged = TRUE,
    hessian_pd = FALSE,
    estimate = 0.5,
    se = 0.1,
    alpha = 0.05
  )
  expect_equal(bad_hess$status, "hessian_not_pd")
  expect_equal(bad_hess$hit_p, NA)

  bad_se <- beezdiscounting:::.dd_power_rep_row(
    converged = TRUE,
    hessian_pd = TRUE,
    estimate = 0.5,
    se = NA_real_,
    alpha = 0.05
  )
  expect_equal(bad_se$status, "se_unusable")
  expect_equal(bad_se$hit_p, NA)

  null_est <- beezdiscounting:::.dd_power_rep_row(
    converged = TRUE,
    hessian_pd = TRUE,
    estimate = 0.001,
    se = 0.5,
    alpha = 0.05
  )
  expect_equal(null_est$status, "ok")
  expect_equal(null_est$hit_p, FALSE)
  expect_equal(null_est$hit_ci, FALSE)
})

test_that(".dd_power_rep_row refers the Wald statistic to a t distribution", {
  # |z| = 2.1: rejects under z (crit 1.96) but not under t with 5 df
  # (crit 2.571); p-value and CI verdicts must agree in both cases.
  z_row <- beezdiscounting:::.dd_power_rep_row(
    converged = TRUE,
    hessian_pd = TRUE,
    estimate = 0.21,
    se = 0.1,
    alpha = 0.05,
    df = Inf
  )
  expect_equal(z_row$hit_p, TRUE)
  expect_equal(z_row$hit_ci, TRUE)

  t_row <- beezdiscounting:::.dd_power_rep_row(
    converged = TRUE,
    hessian_pd = TRUE,
    estimate = 0.21,
    se = 0.1,
    alpha = 0.05,
    df = 5
  )
  expect_equal(t_row$hit_p, FALSE)
  expect_equal(t_row$hit_ci, FALSE)
  expect_equal(t_row$p_value, 2 * stats::pt(-2.1, df = 5))
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

test_that("find_n_discounting validates n_range and n_sim_max", {
  expect_error(
    find_n_discounting(
      target_power = 0.8,
      effect = list(delta_k = 0.5),
      n_range = c(10, 5)
    ),
    "n_range"
  )
  expect_error(
    find_n_discounting(
      target_power = 0.8,
      effect = list(delta_k = 0.5),
      n_sim = 100,
      n_sim_max = 50
    ),
    "n_sim_max"
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

# -----------------------------------------------------------------------------
# Search decision logic under forced disagreement (fake run_batch; no fitting)
# -----------------------------------------------------------------------------

dd_fake_batch <- function(hit_rates_by_n) {
  function(n, batch_size, sim_offset) {
    rate <- hit_rates_by_n[[as.character(n)]]
    if (is.null(rate)) {
      rate <- 0
    }
    hits <- round(batch_size * rate)
    hit <- c(rep(TRUE, hits), rep(FALSE, batch_size - hits))
    tibble::tibble(
      sim = seq_len(batch_size) + sim_offset,
      status = "ok",
      converged = TRUE,
      hessian_pd = TRUE,
      estimate = ifelse(hit, 1, 0),
      se = 0.1,
      statistic = ifelse(hit, 10, 0),
      p_value = ifelse(hit, 1e-6, 1),
      ci_lower = ifelse(hit, 0.8, -0.2),
      ci_upper = ifelse(hit, 1.2, 0.2),
      hit_p = hit,
      hit_ci = hit,
      message = NA_character_
    )
  }
}

test_that("the search confirms a clean minimum", {
  res <- beezdiscounting:::.dd_power_find_n_search(
    dd_fake_batch(list(`4` = 0.1, `7` = 0.2, `8` = 0.99, `12` = 0.99)),
    target_power = 0.8,
    n_range = c(4, 12),
    n_sim = 200,
    n_sim_max = 400,
    verbose = FALSE
  )
  expect_equal(res$n, 8)
  expect_equal(res$status, "confirmed")
  expect_equal(res$uncertain, FALSE)
  expect_true(all(c("n_used", "usable_fraction") %in% names(res$evaluations)))
})

test_that("the lower bound is reconfirmed with fresh replicates before at_lower_bound", {
  res <- beezdiscounting:::.dd_power_find_n_search(
    dd_fake_batch(list(`4` = 0.99, `12` = 0.99)),
    target_power = 0.8,
    n_range = c(4, 12),
    n_sim = 200,
    n_sim_max = 400,
    verbose = FALSE
  )
  expect_equal(res$n, 4)
  expect_equal(res$status, "at_lower_bound")
  expect_equal(res$uncertain, FALSE)
  expect_equal(sum(res$evaluations$n_subjects == 4), 2L)
})

test_that("a lower bound that fails reconfirmation is bisected past, not reported", {
  # First look at 4 reads above, the fresh look reads below: 4 is not
  # reliably above, so [4, 12] is a valid bracket and the search continues
  # upward; the contradictory first look demotes the status to "uncertain".
  calls <- new.env()
  calls$n4 <- 0L
  stateful <- function(n, batch_size, sim_offset) {
    if (n == 4) {
      calls$n4 <- calls$n4 + 1L
      rate <- if (calls$n4 > 1L) 0.1 else 0.99
      return(dd_fake_batch(setNames(list(rate), "4"))(n, batch_size, sim_offset))
    }
    dd_fake_batch(list(`12` = 0.99))(n, batch_size, sim_offset)
  }
  res <- suppressWarnings(beezdiscounting:::.dd_power_find_n_search(
    stateful,
    target_power = 0.8,
    n_range = c(4, 12),
    n_sim = 200,
    n_sim_max = 400,
    verbose = FALSE
  ))
  expect_equal(res$n, 12)
  expect_equal(res$status, "uncertain")
  expect_equal(res$uncertain, TRUE)
  expect_equal(sum(res$evaluations$n_subjects == 4), 2L)
})

test_that("a failed confirmation returns NA with status unresolved", {
  calls <- new.env()
  calls$n8 <- 0L
  stateful <- function(n, batch_size, sim_offset) {
    base <- dd_fake_batch(list(`4` = 0.1, `7` = 0.1, `8` = 0.99, `12` = 0.99))
    if (n == 8) {
      calls$n8 <- calls$n8 + 1L
      if (calls$n8 > 1L) {
        return(dd_fake_batch(list(`8` = 0.1))(n, batch_size, sim_offset))
      }
    }
    base(n, batch_size, sim_offset)
  }
  res <- suppressWarnings(beezdiscounting:::.dd_power_find_n_search(
    stateful,
    target_power = 0.8,
    n_range = c(4, 12),
    n_sim = 200,
    n_sim_max = 400,
    verbose = FALSE
  ))
  expect_equal(res$n, NA_integer_)
  expect_equal(res$status, "unresolved")
  expect_equal(res$uncertain, TRUE)
})
