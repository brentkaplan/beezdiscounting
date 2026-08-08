test_that("lookup21 ships the canonical Kirby & Marakovic (1996) design", {
  lk <- beezdiscounting:::lookup21
  expect_equal(nrow(lk), 21L)
  expect_named(
    lk,
    c(
      "questionid",
      "magnitude",
      "kindiff",
      "k_rank",
      "ss_amount",
      "ll_amount",
      "delay"
    )
  )
  # ladder order is load-bearing (scorer order = Table 1 order)
  expect_equal(
    lk$questionid,
    c(4, 15, 7, 20, 9, 12, 8, 16, 14, 10, 3, 18, 11, 2, 19, 21, 6, 17, 5, 13, 1)
  )
  expect_equal(
    lk$kindiff,
    c(
      0.0007,
      0.0007,
      0.0007,
      0.0032,
      0.0032,
      0.0031,
      0.0057,
      0.0055,
      0.0083,
      0.0089,
      0.0077,
      0.0086,
      0.0160,
      0.0150,
      0.0159,
      0.0250,
      0.0359,
      0.0375,
      0.1333,
      0.1292,
      0.1310
    )
  )
  expect_equal(
    lk$k_rank,
    c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7)
  )
  # per-magnitude structure: 7 each, strictly ascending kindiff
  for (m in c("S", "M", "L")) {
    ks <- lk$kindiff[lk$magnitude == m]
    expect_length(ks, 7L)
    expect_true(all(diff(ks) > 0))
  }
  # Kirby Eq. 1: published k reproduces from amounts/delay within rounding
  k_eq1 <- (lk$ll_amount / lk$ss_amount - 1) / lk$delay
  expect_true(max(abs(k_eq1 - lk$kindiff) / lk$kindiff) < 0.05)
  # spot-check three items verbatim from Table 1
  expect_equal(
    unlist(
      lk[lk$questionid == 4, c("ss_amount", "ll_amount", "delay")],
      use.names = FALSE
    ),
    c(34, 35, 43)
  )
  expect_equal(
    unlist(
      lk[lk$questionid == 1, c("ss_amount", "ll_amount", "delay")],
      use.names = FALSE
    ),
    c(30, 85, 14)
  )
  expect_equal(
    unlist(
      lk[lk$questionid == 10, c("ss_amount", "ll_amount", "delay")],
      use.names = FALSE
    ),
    c(40, 65, 70)
  )
})

test_that("lookup (27-item) is untouched", {
  lk <- beezdiscounting:::lookup
  expect_equal(nrow(lk), 27L)
  expect_equal(lk$questionid[1:4], c(13, 1, 9, 20))
})
