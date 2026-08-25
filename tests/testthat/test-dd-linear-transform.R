# tests/testthat/test-dd-linear-transform.R
describe(".dd_lin_transform()", {
  it("implements Eq. 6: y = ln(1/D - 1) - ln(t) on interior points", {
    d <- c(0.8, 0.5, 0.2); t <- c(1, 7, 30)
    out <- .dd_lin_transform(d, t)
    expect_equal(out$y_lin, log(1 / d - 1) - log(t))
    expect_equal(out$d_used, d)
    expect_equal(out$n_boundary, 0L)
    expect_equal(out$n_clamped, 0L)
    expect_equal(out$log_jac, -log(d - d^2))
  })
  it("clamps D in {0,1} to [eps, 1-eps] and counts them (default)", {
    d <- c(1, 0.5, 0); t <- c(1, 7, 30)
    out <- .dd_lin_transform(d, t, eps = 0.005)
    expect_equal(out$d_used, c(0.995, 0.5, 0.005))
    expect_equal(out$n_boundary, 2L)
    expect_equal(out$n_clamped, 2L)
    expect_equal(out$y_lin, log(1 / c(0.995, 0.5, 0.005) - 1) - log(t))
  })
  it("clamp also pulls interior values inside (eps, 1-eps) and counts those", {
    out <- .dd_lin_transform(c(0.999, 0.5), c(1, 7), eps = 0.005)
    expect_equal(out$d_used, c(0.995, 0.5))
    expect_equal(out$n_boundary, 0L)
    expect_equal(out$n_clamped, 1L)
  })
  it("drop mode returns NA for boundary points and counts drops", {
    out <- .dd_lin_transform(c(1, 0.5, 0), c(1, 7, 30), boundary = "drop")
    expect_equal(is.na(out$y_lin), c(TRUE, FALSE, TRUE))
    expect_equal(out$n_dropped, 2L)
    expect_equal(out$n_clamped, 0L)
  })
  it("error mode aborts on boundary points, names the count", {
    expect_error(.dd_lin_transform(c(1, 0.5), c(1, 7), boundary = "error"), "1 indifference point")
  })
  it("rejects t <= 0 and D outside [0, 1]", {
    expect_error(.dd_lin_transform(0.5, 0), "positive")
    expect_error(.dd_lin_transform(1.2, 1), "\\[0, 1\\]")
  })
})
