# tests/testthat/test-build-fixed-rhs.R

describe("build_fixed_rhs()", {
  it("returns intercept-only ~ 1 with no factors or covariates", {
    f <- build_fixed_rhs()
    expect_s3_class(f, "formula")
    expect_equal(deparse(f[[2]]), "1")
  })

  it("emits a single factor as a main effect", {
    dat <- data.frame(grp = factor(c("a", "b")))
    f <- build_fixed_rhs(factors = "grp", data = dat)
    expect_equal(deparse(f[[2]]), "grp")
  })

  it("emits two factors as additive main effects by default", {
    dat <- data.frame(a = factor(c("x", "y")), b = factor(c("p", "q")))
    f <- build_fixed_rhs(factors = c("a", "b"), data = dat)
    expect_equal(deparse(f[[2]]), "a + b")
  })

  it("includes the interaction when factor_interaction = TRUE", {
    dat <- data.frame(a = factor(c("x", "y")), b = factor(c("p", "q")))
    f <- build_fixed_rhs(factors = c("a", "b"),
      factor_interaction = TRUE, data = dat)
    expect_equal(deparse(f[[2]]), "a * b")
  })

  it("appends continuous covariates as main effects", {
    f <- build_fixed_rhs(continuous_covariates = c("age", "ses"))
    expect_equal(deparse(f[[2]]), "age + ses")
  })

  it("combines a factor and a covariate", {
    dat <- data.frame(grp = factor(c("a", "b")), age = c(20, 30))
    f <- build_fixed_rhs(factors = "grp", continuous_covariates = "age",
      data = dat)
    expect_equal(deparse(f[[2]]), "grp + age")
  })

  it("drops a single-level factor (no contrasts) with a note", {
    dat <- data.frame(grp = factor(c("a", "a")), age = c(20, 30))
    expect_message(
      f <- build_fixed_rhs(factors = "grp", continuous_covariates = "age",
        data = dat),
      regexp = "1 level|removed"
    )
    expect_equal(deparse(f[[2]]), "age")
  })

  it("returns a one-sided formula (no LHS)", {
    f <- build_fixed_rhs()
    expect_true(length(f) == 2L)  # one-sided: c("~", rhs)
  })
})

describe("%||%", {
  it("returns the LHS when not NULL, else the RHS", {
    expect_equal(3L %||% 7L, 3L)
    expect_equal(NULL %||% 7L, 7L)
  })
})
