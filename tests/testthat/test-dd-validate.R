# tests/testthat/test-dd-validate.R

make_ip <- function(y, x = c(1, 7, 30, 90, 180), id = "P1") {
  data.frame(id = id, x = x, y = y, stringsAsFactors = FALSE)
}

describe(".dd_validate_ip", {
  it("passes proportion data [0,1] through unchanged with no warning", {
    dat <- make_ip(c(0.9, 0.6, 0.3, 0.1, 0.0))
    res <- expect_silent(.dd_validate_ip(dat, "y", "x", "id"))
    expect_named(res, c("data", "coercion_info"))
    expect_equal(res$data$y, dat$y)
    expect_identical(res$coercion_info$scale_detected, "proportion")
    expect_equal(res$coercion_info$divided_by, 1)
    expect_equal(res$coercion_info$n_clamped_hi, 0L)
    expect_equal(res$coercion_info$n_clamped_lo, 0L)
  })

  it("detects percent (max > 1.5), divides by 100, and WARNS", {
    dat <- make_ip(c(90, 60, 30, 10, 0))
    expect_warning(
      res <- .dd_validate_ip(dat, "y", "x", "id"),
      regexp = "percent|100"
    )
    expect_equal(res$data$y, c(0.9, 0.6, 0.3, 0.1, 0))
    expect_identical(res$coercion_info$scale_detected, "percent")
    expect_equal(res$coercion_info$divided_by, 100)
  })

  it("divides by a supplied larger-later reward (amount) and WARNS", {
    dat <- make_ip(c(900, 600, 300, 100, 0))
    expect_warning(
      res <- .dd_validate_ip(dat, "y", "x", "id", ll = 1000,
        response_scale = "amount"),
      regexp = "amount|1000|larger-later"
    )
    expect_equal(res$data$y, c(0.9, 0.6, 0.3, 0.1, 0))
    expect_identical(res$coercion_info$scale_detected, "amount")
    expect_equal(res$coercion_info$divided_by, 1000)
  })

  it("clamps mild out-of-range after coercion and WARNS with named counts", {
    # proportion scale but two strays: 1.02 -> 1 (hi), -0.01 -> 0 (lo)
    dat <- make_ip(c(1.02, 0.6, -0.01, 0.1, 0.0))
    expect_warning(
      res <- .dd_validate_ip(dat, "y", "x", "id"),
      regexp = "[Cc]lamped"
    )
    expect_equal(res$data$y, c(1, 0.6, 0, 0.1, 0))
    expect_equal(res$coercion_info$n_clamped_hi, 1L)
    expect_equal(res$coercion_info$n_clamped_lo, 1L)
  })

  it("errors when id/x/y columns are missing", {
    bad <- data.frame(subject = "P1", delay = 1, ip = 0.5)
    expect_error(.dd_validate_ip(bad, "y", "x", "id"), regexp = "not found|missing")
  })

  it("errors when y has NA after coercion", {
    dat <- make_ip(c(0.9, NA, 0.3, 0.1, 0.0))
    expect_error(.dd_validate_ip(dat, "y", "x", "id"), regexp = "NA")
  })

  it("remaps caller column names via y_var/x_var/id_var", {
    dat <- data.frame(
      PID = "P1", delay = c(1, 7, 30), prop = c(0.9, 0.5, 0.1),
      stringsAsFactors = FALSE
    )
    res <- .dd_validate_ip(dat, y_var = "prop", x_var = "delay", id_var = "PID")
    expect_named(res$data, c("id", "x", "y"))
    expect_equal(res$data$y, c(0.9, 0.5, 0.1))
    expect_equal(res$data$x, c(1, 7, 30))
  })

  it("retains extra factor/covariate columns alongside id/x/y (row-coherence)", {
    dat <- data.frame(
      id = "P1",
      x = c(1, 7, 30),
      y = c(0.9, 0.5, 0.1),
      condition = factor(c("A", "A", "B"), levels = c("A", "B")),
      age = c(20, 20, 20),
      stringsAsFactors = FALSE
    )
    res <- .dd_validate_ip(dat, "y", "x", "id",
      extra_cols = c("condition", "age"))
    # id/x/y come first, then the retained columns in requested order
    expect_named(res$data, c("id", "x", "y", "condition", "age"))
    # factor class + levels preserved verbatim (not coerced to character)
    expect_s3_class(res$data$condition, "factor")
    expect_identical(levels(res$data$condition), c("A", "B"))
    expect_equal(res$data$age, c(20, 20, 20))
  })

  it("errors when a requested extra column is missing", {
    dat <- make_ip(c(0.9, 0.6, 0.3, 0.1, 0.0))
    expect_error(
      .dd_validate_ip(dat, "y", "x", "id", extra_cols = "nope"),
      regexp = "not found"
    )
  })

  it("passes proportion data with max == 1.0 through silently (not mis-detected as percent)", {
    dat <- make_ip(c(1.0, 0.8, 0.5, 0.2, 0.0))
    res <- expect_silent(.dd_validate_ip(dat, "y", "x", "id"))
    expect_identical(res$coercion_info$scale_detected, "proportion")
    expect_equal(res$coercion_info$divided_by, 1)
    expect_equal(res$data$y, c(1.0, 0.8, 0.5, 0.2, 0.0))
  })
})
