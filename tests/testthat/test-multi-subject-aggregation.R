# Regression tests: check_unsystematic(), calc_aucs(), and prop_ss() must handle
# multi-subject input per subject rather than computing one result and recycling it
# across unique(id). See R/utils.R and R/mcq.R.

test_that("check_unsystematic returns a distinct verdict per subject", {
  d <- rbind(
    data.frame(
      id = "a",
      x = c(1, 7, 30, 90, 180, 365),
      y = c(0.95, 0.81, 0.55, 0.32, 0.17, 0.08)
    ), # orderly -> passes
    data.frame(
      id = "b",
      x = c(1, 7, 30, 90, 180, 365),
      y = c(0.20, 0.90, 0.10, 0.80, 0.10, 0.70)
    ) # bouncy -> fails C1
  )
  out <- check_unsystematic(d)

  expect_equal(nrow(out), 2L)
  expect_setequal(out$id, c("a", "b"))
  expect_true(out$c1_pass[out$id == "a"])
  expect_false(out$c1_pass[out$id == "b"])
})

test_that("check_unsystematic single-subject result is unchanged", {
  d <- data.frame(id = "a", x = c(1, 7, 30), y = c(0.9, 0.5, 0.1))
  out <- check_unsystematic(d)

  expect_equal(nrow(out), 1L)
  expect_true(out$c1_pass)
  expect_true(out$c2_pass)
})

test_that("calc_aucs computes a distinct AUC per subject (no recycling)", {
  d <- rbind(
    data.frame(id = "a", x = c(1, 10, 100), y = c(0.9, 0.5, 0.1)),
    data.frame(id = "b", x = c(1, 10, 100), y = c(0.5, 0.3, 0.2))
  )
  out <- calc_aucs(d)

  expect_equal(nrow(out), 2L)
  expect_setequal(out$id, c("a", "b"))
  expect_false(isTRUE(all.equal(
    out$auc_regular[out$id == "a"],
    out$auc_regular[out$id == "b"]
  )))

  # per-subject value matches calc_aucs() on that subject alone
  solo_a <- calc_aucs(d[d$id == "a", ])
  expect_equal(out$auc_regular[out$id == "a"], solo_a$auc_regular)
})

test_that("prop_ss returns proportions within [0, 1] for multiple respondents", {
  set.seed(42)
  many <- generate_data_mcq(n_ids = 10) # old sum/3 denominator could reach ~10
  out <- prop_ss(many)

  expect_equal(nrow(out), 9L) # nine k-ranks
  expect_true(all(out$prop_ss >= 0 & out$prop_ss <= 1))
})

test_that("calc_aucs accepts id-less data as a single subject", {
  d <- data.frame(x = c(1, 7, 30, 90), y = c(0.9, 0.6, 0.3, 0.1))
  out <- calc_aucs(d)

  expect_equal(nrow(out), 1L)
  expect_false("id" %in% names(out))
  expect_true(all(c("auc_regular", "auc_log10", "auc_ord") %in% names(out)))
  expect_true(is.finite(out$auc_regular))
})

test_that("results_dd(method = 'mean') still returns AUC columns", {
  d <- rbind(
    data.frame(id = "a", x = c(1, 7, 30, 90), y = c(0.90, 0.60, 0.30, 0.10)),
    data.frame(id = "b", x = c(1, 7, 30, 90), y = c(0.80, 0.50, 0.20, 0.05))
  )
  res <- results_dd(fit_dd(d, equation = "mazur", method = "mean"))

  expect_true(all(c("auc_regular", "auc_log10", "auc_ord") %in% names(res)))
  expect_true(is.finite(res$auc_regular[[1]]))
})

test_that("prop_ss pools across respondents instead of dropping all but one", {
  # subject a always chooses SS (0), subject b always LL (1); pooled SS share = 0.5.
  # The old match() reduced this to one respondent, returning 1.0 at every rank.
  two <- rbind(
    data.frame(subjectid = "a", questionid = 1:27, response = 0L),
    data.frame(subjectid = "b", questionid = 1:27, response = 1L)
  )
  out <- prop_ss(two)

  expect_equal(nrow(out), 9L)
  expect_true(all(out$prop_ss == 0.5))
})

test_that("check_unsystematic ignores rows with a missing id", {
  # Subject "a" is perfectly orderly; NA-id rows must not leak into its subset
  # (base subsetting on `id == i` would otherwise inject synthetic NA rows).
  solo <- data.frame(id = "a", x = c(1, 30, 365), y = c(0.9, 0.5, 0.1))
  d <- rbind(
    solo,
    data.frame(id = NA_character_, x = c(7, 90), y = c(0.4, 0.3))
  )
  out <- check_unsystematic(d)

  expect_equal(nrow(out), 1L)
  expect_equal(out$id, "a")
  expect_true(out$c1_pass)
  expect_true(out$c2_pass)
})

test_that("calc_aucs ignores rows with a missing id", {
  solo <- data.frame(id = "a", x = c(1, 10, 100), y = c(0.9, 0.5, 0.1))
  d <- rbind(
    solo,
    data.frame(id = NA_character_, x = c(5, 50), y = c(0.6, 0.2))
  )
  out <- calc_aucs(d)

  expect_equal(nrow(out), 1L)
  expect_equal(out$id, "a")
  expect_equal(out$auc_regular, calc_aucs(solo)$auc_regular)
})
