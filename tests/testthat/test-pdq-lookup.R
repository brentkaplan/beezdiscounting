test_that("lookup_pdq ships the canonical Madden et al. (2009) design", {
  lk <- beezdiscounting:::lookup_pdq
  expect_equal(nrow(lk), 30L)
  expect_named(
    lk,
    c(
      "questionid",
      "block",
      "h_rank",
      "overall_rank",
      "sc_amount",
      "lu_amount",
      "prob",
      "theta",
      "hindiff"
    )
  )
  expect_identical(lk$questionid, 1:30)
  expect_identical(lk$block, rep(1:3, each = 10L))
  expect_identical(lk$h_rank, rep(1:10, times = 3L))
  # pooled-ladder order: exact-rational h ascending, ties broken by
  # questionid (q1&q21, q4&q24, q5&q15, q16&q26 are exact ties)
  expect_identical(
    lk$questionid[order(lk$overall_rank)],
    c(
      11L,
      1L,
      21L,
      12L,
      22L,
      2L,
      23L,
      13L,
      3L,
      14L,
      4L,
      24L,
      5L,
      15L,
      25L,
      6L,
      16L,
      26L,
      7L,
      17L,
      27L,
      28L,
      18L,
      8L,
      9L,
      19L,
      29L,
      10L,
      20L,
      30L
    )
  )
  # hindiff is (weakly) nondecreasing along the pooled order up to
  # floating-point noise at the exact-tie pairs
  h_pooled <- lk$hindiff[order(lk$overall_rank)]
  expect_true(all(diff(h_pooled) > -1e-9))
  expect_equal(lk$sc_amount, rep(c(20, 40, 40), each = 10))
  expect_equal(lk$lu_amount, rep(c(80, 100, 60), each = 10))
  expect_equal(
    lk$prob,
    c(
      .10,
      .13,
      .17,
      .20,
      .25,
      .33,
      .50,
      .67,
      .75,
      .83,
      .18,
      .22,
      .29,
      .33,
      .40,
      .50,
      .67,
      .80,
      .86,
      .91,
      .40,
      .46,
      .55,
      .60,
      .67,
      .75,
      .86,
      .92,
      .95,
      .97
    )
  )
  # theta and hindiff are exact functions of the printed design values
  expect_equal(lk$theta, (1 - lk$prob) / lk$prob, tolerance = 1e-12)
  expect_equal(
    lk$hindiff,
    (lk$lu_amount / lk$sc_amount - 1) / lk$theta,
    tolerance = 1e-12
  )
  # h strictly ascending within each block (ladder order is load-bearing)
  for (b in 1:3) {
    hs <- lk$hindiff[lk$block == b]
    expect_length(hs, 10L)
    expect_true(all(diff(hs) > 0))
  }
  # printed Table 2 values reproduce to ~1% relative error. NOT a
  # rounding-equality check: printed items 13, 17, and 25 are last-digit
  # inconsistent with the computed h (published-table artifacts); the
  # Gray lookup files are the precision oracle (see test-pdq-score.R).
  printed <- c(
    0.33,
    0.45,
    0.61,
    0.75,
    1,
    1.48,
    3,
    6.09,
    9,
    14.65,
    0.33,
    0.42,
    0.62,
    0.74,
    1,
    1.5,
    3.04,
    6,
    9.21,
    15.17,
    0.33,
    0.43,
    0.61,
    0.75,
    1.01,
    1.5,
    3.07,
    5.75,
    9.5,
    16.17
  )
  expect_true(max(abs(lk$hindiff - printed) / printed) < 0.02)
  # spot-check the ladder maxima to lookup-table precision
  expect_equal(
    lk$hindiff[c(10, 20, 30)],
    c(14.6470588235, 15.1666666667, 16.1666666667),
    tolerance = 1e-9
  )
})

test_that("lookup (27) and lookup21 are untouched", {
  expect_equal(nrow(beezdiscounting:::lookup), 27L)
  expect_equal(nrow(beezdiscounting:::lookup21), 21L)
  expect_equal(beezdiscounting:::lookup21$questionid[1:4], c(4, 15, 7, 20))
})
