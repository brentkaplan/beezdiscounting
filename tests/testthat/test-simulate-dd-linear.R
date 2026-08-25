describe("simulate_dd_linear()", {
  it("returns long data in package format with D strictly inside (0,1)", {
    s <- simulate_dd_linear(n_subjects = 5, delays = c(7, 30), mu = -5, sigma2 = 1, g = 4, seed = 1)
    expect_named(s, c("id", "condition", "x", "y"))
    expect_equal(nrow(s), 10L); expect_true(all(s$y > 0 & s$y < 1))
    expect_equal(levels(s$condition), "C1")
  })
  it("D = 1/(1 + exp(y) t) inverts the transform exactly", {
    s <- simulate_dd_linear(3, c(7, 30, 90), mu = -5, sigma2 = 1, g = 4, seed = 2, attach_truth = TRUE)
    tr <- .dd_lin_transform(s$y, s$x)
    expect_equal(1 / (1 + exp(tr$y_lin) * s$x), s$y)
  })
  it("named mu gives named conditions; n_subjects recycles per condition", {
    s <- simulate_dd_linear(c(2, 3), c(7, 30), mu = c(EFT = -7, NCC = -5.5), sigma2 = 2, g = 10, seed = 3)
    expect_equal(as.integer(table(s$condition)), c(4L, 6L)); expect_equal(levels(s$condition), c("EFT", "NCC"))
  })
  it("seed makes it reproducible and restores the global RNG state", {
    set.seed(99); invisible(runif(1))
    state_before <- .Random.seed
    a <- simulate_dd_linear(4, c(7, 30), -5, 1, 4, seed = 7)
    expect_identical(.Random.seed, state_before)          # state untouched by the seeded call
    b <- simulate_dd_linear(4, c(7, 30), -5, 1, 4, seed = 7)
    expect_identical(a, b)
  })
  it("RECOVERY: fit_dd_linear recovers mu, sigma2, g at paper Sec 4.3 values", {
    skip_on_cran()
    s <- simulate_dd_linear(n_subjects = 200, delays = c(30, 90, 180, 365, 1095, 1825, 3650),
                            mu = c(-7, -5.5), sigma2 = 2, g = 10.4, seed = 42)
    f <- fit_dd_linear(s, factors = "condition")
    expect_equal(unname(f$re$mu), c(-7, -5.5), tolerance = 0.1)
    expect_equal(f$re$sigma2, 2, tolerance = 0.08)
    expect_equal(f$re$g, 10.4, tolerance = 0.15)
  })
  it("TYPE I: all-equal F rejects 5% under H0 (2000 predeclared reps, 99% binomial band)", {
    skip_on_cran()
    # Predeclared: seeds 100001..102000; band = 0.05 +/- 2.576 * sqrt(0.05 * 0.95 / 2000) = [0.037, 0.063].
    # A failure is a bug to investigate (check df2 and that SSR_{Z|X_full}, not SSE_Z, is the
    # denominator) -- never re-seed or widen the band.
    rej <- vapply(seq_len(2000), function(i) {
      s <- simulate_dd_linear(35, c(30, 90, 180, 365, 1095, 1825, 3650), mu = c(-5.5, -5.5),
                              sigma2 = 2, g = 10.4, seed = 100000 + i)
      # .dd_lin_ftest() is what anova() wraps (Task 6); kept internal here so this file is self-contained.
      .dd_lin_ftest(fit_dd_linear(s, factors = "condition")$re)$p_value < 0.05
    }, logical(1))
    expect_gte(mean(rej), 0.037); expect_lte(mean(rej), 0.063)
  })
})
