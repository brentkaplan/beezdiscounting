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

  it("normalizes k + phi ~ 1 to a 2-D block; covariance_structure sets the class", {
    res <- .dd_normalize_re(k + phi ~ 1, covariance_structure = "pdSymm")
    expect_length(res$blocks, 1L)
    b <- res$blocks[[1]]
    expect_identical(b$param, c("k", "phi"))
    expect_identical(b$terms, "(Intercept)")
    expect_identical(b$pdmat_class, "pdSymm")
    expect_identical(b$dim, 2L)
  })

  it("canonicalizes phi + k ~ 1 to (k, phi) order", {
    res <- .dd_normalize_re(phi + k ~ 1, covariance_structure = "pdDiag")
    expect_identical(res$blocks[[1]]$param, c("k", "phi"))
    expect_identical(res$blocks[[1]]$pdmat_class, "pdDiag")
  })

  it("defaults covariance_structure to pdSymm for the 2-D block", {
    res <- .dd_normalize_re(k + phi ~ 1)
    expect_identical(res$blocks[[1]]$pdmat_class, "pdSymm")
  })

  it("rejects an unknown covariance_structure", {
    expect_error(.dd_normalize_re(k + phi ~ 1, covariance_structure = "pdLogChol"),
      regexp = "pdDiag|pdSymm")
  })

  it("rejects random slopes (k ~ 1 + condition)", {
    expect_error(.dd_normalize_re(k ~ 1 + condition),
      regexp = "intercept-only|slope|single")
  })

  it("rejects a transformed LHS; bare symbols required (B10)", {
    expect_error(.dd_normalize_re(log(k) ~ 1), "bare symbol")
    expect_error(.dd_normalize_re(I(k) ~ 1), "bare symbol")
    expect_silent(.dd_normalize_re(k ~ 1))
    expect_silent(.dd_normalize_re(`k` ~ 1))
  })

  it("rejects a phi-only LHS (phi ~ 1) — k is required", {
    expect_error(.dd_normalize_re(phi ~ 1), regexp = "must include.*k|include the .*k")
  })

  it("rejects a third LHS parameter (k + phi + s ~ 1)", {
    expect_error(.dd_normalize_re(k + phi + s ~ 1), regexp = "k.*phi|only.*k.*phi")
  })

  it("rejects a duplicated LHS parameter (k + k ~ 1)", {
    expect_error(.dd_normalize_re(k + k ~ 1), regexp = "duplicat|repeat")
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
