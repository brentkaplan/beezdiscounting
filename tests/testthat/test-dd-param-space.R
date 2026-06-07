# tests/testthat/test-dd-param-space.R

describe(".dd_param_registry", {
  it("has entries keyed by k, phi, s with the required fields", {
    expect_setequal(names(.dd_param_registry), c("k", "phi", "s"))
    for (key in names(.dd_param_registry)) {
      reg <- .dd_param_registry[[key]]
      expect_true(all(c(
        "canonical", "description", "constraint",
        "valid_scales", "default_scale"
      ) %in% names(reg)))
      expect_identical(reg$canonical, key)
      expect_true(reg$default_scale %in% reg$valid_scales)
    }
  })

  it("declares k on natural/log/log10 defaulting to natural", {
    expect_setequal(.dd_param_registry$k$valid_scales,
      c("natural", "log", "log10"))
    expect_identical(.dd_param_registry$k$default_scale, "natural")
  })
})

describe(".dd_transform_est_se", {
  it("is a no-op when from == to or to == 'internal'", {
    r <- .dd_transform_est_se(2, 0.5, "natural", "natural")
    expect_equal(r$estimate, 2)
    expect_equal(r$se, 0.5)
    r2 <- .dd_transform_est_se(2, 0.5, "log", "internal")
    expect_equal(r2$estimate, 2)
    expect_equal(r2$se, 0.5)
  })

  it("round-trips natural <-> log10 (delta method on SE)", {
    fwd <- .dd_transform_est_se(0.01, 0.002, "natural", "log10")
    expect_equal(fwd$estimate, log10(0.01), tolerance = 1e-12)
    expect_equal(fwd$se, 0.002 / (0.01 * log(10)), tolerance = 1e-12)
    back <- .dd_transform_est_se(fwd$estimate, fwd$se, "log10", "natural")
    expect_equal(back$estimate, 0.01, tolerance = 1e-10)
    expect_equal(back$se, 0.002, tolerance = 1e-10)
  })

  it("round-trips natural <-> log (delta method on SE)", {
    fwd <- .dd_transform_est_se(0.01, 0.002, "natural", "log")
    expect_equal(fwd$estimate, log(0.01), tolerance = 1e-12)
    expect_equal(fwd$se, 0.002 / 0.01, tolerance = 1e-12)
    back <- .dd_transform_est_se(fwd$estimate, fwd$se, "log", "natural")
    expect_equal(back$estimate, 0.01, tolerance = 1e-10)
    expect_equal(back$se, 0.002, tolerance = 1e-10)
  })

  it("converts log <-> log10 by ln(10)", {
    r <- .dd_transform_est_se(log(5), 0.3, "log", "log10")
    expect_equal(r$estimate, log(5) / log(10), tolerance = 1e-12)
    expect_equal(r$se, 0.3 / log(10), tolerance = 1e-12)
  })

  it("returns NA estimate/se for non-positive values into log/log10", {
    r <- .dd_transform_est_se(-1, 0.2, "natural", "log10")
    expect_true(is.na(r$estimate))
    expect_true(is.na(r$se))
  })

  it("errors on an unsupported transform pair", {
    expect_error(.dd_transform_est_se(1, 0.1, "natural", "bogus"))
  })
})

describe(".dd_transform_coef_table", {
  it("transforms only core k/phi/s rows and leaves others untouched", {
    tbl <- data.frame(
      term = c("k", "phi", "s", "(Intercept)_z"),
      estimate = c(0.01, 8, 1, 99),
      std.error = c(0.002, 1, 0.1, 1),
      stringsAsFactors = FALSE
    )
    out <- .dd_transform_coef_table(tbl, report_space = "log10",
      internal_space = "natural")
    # core rows -> log10
    expect_equal(out$estimate[out$term == "k"], log10(0.01), tolerance = 1e-12)
    expect_equal(out$estimate[out$term == "phi"], log10(8), tolerance = 1e-12)
    # non-core row passes through unchanged
    expect_equal(out$estimate[out$term == "(Intercept)_z"], 99)
    # internal estimate preserved for core rows
    expect_equal(out$estimate_internal[out$term == "k"], 0.01)
  })

  it("matches the retargeted regex for prefixed terms (k_log, phi_log)", {
    tbl <- data.frame(
      term = c("k_log", "phi_log", "kappa_not_core"),
      estimate = c(0.02, 10, 5),
      std.error = c(0.001, 1, 1),
      stringsAsFactors = FALSE
    )
    out <- .dd_transform_coef_table(tbl, report_space = "log",
      internal_space = "natural")
    expect_equal(out$estimate[out$term == "k_log"], log(0.02), tolerance = 1e-12)
    expect_equal(out$estimate[out$term == "phi_log"], log(10), tolerance = 1e-12)
    # 'kappa_not_core' must NOT match ^k($|_|:) (next char is 'a', not $, _, or :)
    expect_equal(out$estimate[out$term == "kappa_not_core"], 5)
  })

  it("matches the renamed display terms with a ':' on the k branch", {
    tbl <- data.frame(
      term = c("k:(Intercept)", "k:condition_B", "kappa:not_core"),
      estimate = c(0.01, 0.02, 5),
      std.error = c(0.002, 0.001, 1),
      stringsAsFactors = FALSE
    )
    out <- .dd_transform_coef_table(tbl, report_space = "log10",
      internal_space = "natural")
    # both renamed k display terms are treated as core k-rows
    expect_equal(out$estimate[out$term == "k:(Intercept)"], log10(0.01),
      tolerance = 1e-12)
    expect_equal(out$estimate[out$term == "k:condition_B"], log10(0.02),
      tolerance = 1e-12)
    expect_equal(out$estimate_internal[out$term == "k:(Intercept)"], 0.01)
    # 'kappa:not_core' must NOT match ^k($|_|:) (next char is 'a')
    expect_equal(out$estimate[out$term == "kappa:not_core"], 5)
  })
})

describe("%||%", {
  it("returns x when x is not NULL", {
    expect_equal(3L %||% 7L, 3L)
    expect_equal("a" %||% "b", "a")
  })

  it("returns y when x is NULL", {
    expect_equal(NULL %||% 7, 7)
    expect_equal(NULL %||% "fallback", "fallback")
  })
})
