# Audit 2026-09-06 F-BZ3-1: the Johnson & Bickel (2008) thresholds are
# proportions of the larger-later reward, so they must be applied on the y scale
# (c1 * ll, c2 * ll). Before the fix `ll` cancelled out and amount-scale data
# were judged against the bare proportions.

test_that("check_unsystematic applies C1 on the y scale (amount-scale data)", {
  # LL = 100: a 15-unit bounce is 15% of LL -> passes C1 (J&B criterion 1)
  d <- tibble::tibble(id = 1, x = c(1, 7, 30, 90), y = c(90, 50, 65, 10))
  expect_true(check_unsystematic(d, ll = 100)$c1_pass)
  # a 25-unit bounce is 25% of LL -> fails C1
  d$y[3] <- 75
  expect_false(check_unsystematic(d, ll = 100)$c1_pass)
})

test_that("check_unsystematic applies C2 on the y scale (amount-scale data)", {
  # LL = 100: an 8-unit first-to-last drop is < 10% of LL -> fails C2
  d <- tibble::tibble(id = 1, x = c(1, 30, 90), y = c(90, 85, 82))
  expect_false(check_unsystematic(d, ll = 100)$c2_pass)
  # a 30-unit drop passes
  d$y[3] <- 60
  expect_true(check_unsystematic(d, ll = 100)$c2_pass)
})

test_that("check_unsystematic verdicts are invariant to the unit of y and ll", {
  d <- tibble::tibble(
    id = rep(c("a", "b", "c"), each = 4), x = rep(c(1, 7, 30, 90), 3),
    y = c(0.9, 0.5, 0.65, 0.1, 0.9, 0.85, 0.84, 0.82, 0.95, 0.7, 0.4, 0.2)
  )
  prop <- check_unsystematic(d, ll = 1)
  d_amt <- d
  d_amt$y <- d_amt$y * 250
  amt <- check_unsystematic(d_amt, ll = 250)
  expect_identical(prop$c1_pass, amt$c1_pass)
  expect_identical(prop$c2_pass, amt$c2_pass)
  # and for a tiny unit (the tolerance scales with ll)
  d_small <- d
  d_small$y <- d_small$y * 1e-9
  small <- check_unsystematic(d_small, ll = 1e-9)
  expect_identical(prop$c1_pass, small$c1_pass)
  expect_identical(prop$c2_pass, small$c2_pass)
  expect_false(check_unsystematic(tibble::tibble(id = 1, y = c(1e-9, 1e-9)), ll = 1e-9)$c2_pass)
})

test_that("check_unsystematic thresholds are inclusive at exactly c1 / c2 of ll", {
  # exactly a 10% decline passes C2 (strict: fail only when the drop is < c2 * ll)
  expect_true(check_unsystematic(tibble::tibble(id = 1, y = c(1, 0.9)), ll = 1)$c2_pass)
  expect_true(check_unsystematic(tibble::tibble(id = 1, y = c(100, 90)), ll = 100)$c2_pass)
  # exactly a 20% bounce passes C1
  expect_true(check_unsystematic(tibble::tibble(id = 1, y = c(0.5, 0.7, 0.1)), ll = 1)$c1_pass)
  expect_true(check_unsystematic(tibble::tibble(id = 1, y = c(50, 70, 10)), ll = 100)$c1_pass)
})

test_that("check_unsystematic validates ll, c1 and c2", {
  d <- tibble::tibble(id = 1, y = c(1, 0.5))
  expect_error(check_unsystematic(d, ll = 0), "ll")
  expect_error(check_unsystematic(d, ll = -1), "ll")
  expect_error(check_unsystematic(d, ll = c(1, 2)), "ll")
  expect_error(check_unsystematic(d, ll = NA_real_), "ll")
  expect_error(check_unsystematic(d, c1 = -0.1), "c1")
  expect_error(check_unsystematic(d, c2 = NA_real_), "c2")
})
