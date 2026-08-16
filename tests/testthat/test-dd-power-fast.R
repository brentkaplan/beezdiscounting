# Fast (no-model-fit) tests for power_discounting()/find_n_discounting():
# argument validation, Wilson interval, replicate classification, and the
# mocked bisection search. Split out of the heavy test-dd-power.R
# (BEEZ_FULL_TESTS-gated) so they still run in every ungated test run (CI
# test-coverage, local devtools::test()). See helper-full-tests.R.

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
