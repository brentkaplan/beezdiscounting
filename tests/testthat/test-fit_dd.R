.fit_dd_demo <- function() {
  data.frame(
    id = rep(1:2, each = 6),
    x = rep(c(1, 7, 30, 90, 180, 365), 2),
    y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05,
          0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
  )
}

describe("fit_dd()", {
  it("returns a length-3 fit_dd object for the pooled method", {
    fit <- fit_dd(.fit_dd_demo(), equation = "mazur", method = "pooled")
    expect_s3_class(fit, "fit_dd")
    expect_length(fit, 3L)
    expect_s3_class(fit[[1]], "nls")        # single pooled nlsLM model
    expect_s3_class(fit[[2]]$id, "factor")  # id coerced to factor
    expect_identical(fit[[3]], "pooled")
  })

  it("returns one model for the mean method", {
    fit <- fit_dd(.fit_dd_demo(), equation = "mazur", method = "mean")
    expect_s3_class(fit[[1]], "nls")
    expect_identical(fit[[3]], "mean")
  })

  it("returns one safely() result per id for two stage", {
    fit <- fit_dd(.fit_dd_demo(), equation = "mazur", method = "two stage")
    expect_type(fit[[1]], "list")
    expect_named(fit[[1]], c("1", "2"))
    expect_true(all(c("result", "error") %in% names(fit[[1]][["1"]])))
    expect_s3_class(fit[[1]][["1"]]$result, "nls")
  })

  it("accepts the exponential equation", {
    fit <- fit_dd(.fit_dd_demo(), equation = "exponential", method = "pooled")
    expect_s3_class(fit[[1]], "nls")
  })

  it("accepts mixed-case equation and method", {
    fit <- fit_dd(.fit_dd_demo(), equation = "Hyperbolic", method = "Two Stage")
    expect_identical(fit[[3]], "Two Stage")
  })

  it("rejects an unknown equation", {
    expect_error(fit_dd(.fit_dd_demo(), equation = "weibull", method = "pooled"))
  })

  it("rejects an unknown method", {
    expect_error(fit_dd(.fit_dd_demo(), equation = "mazur", method = "bayes"))
  })
})

describe("results_dd()", {
  it("returns the documented tidy columns for the pooled method", {
    fit <- fit_dd(.fit_dd_demo(), equation = "mazur", method = "pooled")
    out <- results_dd(fit)
    expect_s3_class(out, "tbl_df")
    expect_true(all(c("method", "term", "estimate", "std.error", "statistic",
                      "p.value", "R2", "conf_low", "conf_high") %in% names(out)))
    expect_identical(out$term[1], "k")
    expect_gt(out$estimate[1], 0)
  })

  it("adds AUC columns for the mean method", {
    fit <- fit_dd(.fit_dd_demo(), equation = "mazur", method = "mean")
    out <- results_dd(fit)
    expect_true(any(grepl("^auc", names(out))))
  })

  it("returns one row per id for the two-stage method", {
    fit <- fit_dd(.fit_dd_demo(), equation = "mazur", method = "two stage")
    out <- results_dd(fit)
    expect_setequal(unique(out$id), c("1", "2"))
    expect_true(all(c("method", "id", "estimate", "R2",
                      "conf_low", "conf_high") %in% names(out)))
  })

  it("errors on a non-fit_dd object", {
    expect_error(results_dd(list(1, 2, 3)))
  })
})
