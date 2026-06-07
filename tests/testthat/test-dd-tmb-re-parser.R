# tests/testthat/test-dd-tmb-re-parser.R

describe(".dd_normalize_re", {
  it("normalizes the default k ~ 1 to a single intercept block", {
    res <- .dd_normalize_re(k ~ 1)
    expect_named(res, c("source", "blocks"))
    expect_identical(res$source, "formula")
    expect_length(res$blocks, 1L)
    b <- res$blocks[[1]]
    expect_identical(b$param, "k")
    expect_identical(b$terms, "(Intercept)")
    expect_identical(b$pdmat_class, "pdDiag")
    expect_identical(b$dim, 1L)
  })

  it("rejects random slopes (k ~ 1 + condition)", {
    expect_error(.dd_normalize_re(k ~ 1 + condition),
      regexp = "intercept-only|slope|single")
  })

  it("rejects a non-k LHS (phi ~ 1)", {
    expect_error(.dd_normalize_re(phi ~ 1), regexp = "k")
  })

  it("rejects a two-parameter LHS (k + phi ~ 1)", {
    expect_error(.dd_normalize_re(k + phi ~ 1), regexp = "k")
  })

  it("rejects an intercept-suppressed formula (k ~ 0 + condition)", {
    expect_error(.dd_normalize_re(k ~ 0 + condition),
      regexp = "intercept-only|slope|single")
  })

  it("rejects non-formula input (character vector)", {
    expect_error(.dd_normalize_re(c("k")), regexp = "formula")
  })

  it("rejects a one-sided formula (~ 1)", {
    expect_error(.dd_normalize_re(~ 1), regexp = "two-sided|k ~ 1")
  })
})
