skip_on_cran(); skip_if_not_installed("TMB")
.cdm <- new.env(parent = emptyenv())
.get_desc_fit <- function() {
  if (!exists("fit", envir = .cdm)) {
    dat <- .choice_desc_fixture(n_subjects = 80, seed = 31)
    .cdm$fit <- fit_dd_choice(dat, mode = "descriptive", verbose = 0)
  }
  .cdm$fit
}

describe("beezdiscounting_choice S3 (descriptive)", {
  it("coef/fixef return theta + log_sd_re + cor_re (no structural params)", {
    co <- coef(.get_desc_fit())
    expect_true(all(c("theta", "log_sd_re", "cor_re") %in% names(co)))
    expect_false(any(c("beta_k", "log_gamma") %in% names(co)))
    expect_identical(nlme::fixef(.get_desc_fit()), co)
  })
  it("ranef returns per-subject slopes (id, b_mag, b_delay)", {
    re <- ranef(.get_desc_fit())
    expect_named(re, c("id", "b_mag", "b_delay"))
    expect_equal(nrow(re), .get_desc_fit()$param_info$n_subjects)
  })
  it("VarCorr returns the 2-RE SD/correlation structure", {
    vc <- nlme::VarCorr(.get_desc_fit())
    expect_true(is.data.frame(vc) || is.matrix(vc))
    txt <- paste(capture.output(print(vc)), collapse = " ")
    expect_match(txt, "mag|delay|Corr|SD|StdDev", ignore.case = TRUE)
  })
  it("tidy reports theta rows on the logit (identity) scale, NOT exponentiated", {
    td <- tidy(.get_desc_fit())
    expect_named(td, c("term", "estimate", "std.error", "statistic",
                       "p.value", "component", "estimate_scale", "term_display"))
    mag <- td[grepl("ll_amount|mag", td$term), ]
    expect_true(nrow(mag) >= 1)
    td_nat <- tidy(.get_desc_fit(), report_space = "natural")
    expect_equal(td_nat$estimate[td_nat$component == "fixed"],
                 td$estimate[td$component == "fixed"], tolerance = 1e-10)
  })
  it("predict(type='prob') returns P(LL) in (0,1) of length nobs", {
    pr <- predict(.get_desc_fit(), type = "prob")
    p <- pr$.prob
    expect_length(p, nrow(.get_desc_fit()$data))
    expect_true(all(p > 0 & p < 1))
  })
  it("predict population level (b=0) runs without an id column", {
    nd <- unique(.get_desc_fit()$data[, c("ss_amount", "ll_amount", "delay")])
    pr <- predict(.get_desc_fit(), newdata = nd, type = "prob", level = "population")
    expect_true(all(pr$.prob > 0 & pr$.prob < 1))
  })
  it("predict(level='subject') aborts on unknown ids", {
    nd <- .get_desc_fit()$data[1:3, ]
    nd$id <- "totally-new-subject"
    expect_error(predict(.get_desc_fit(), newdata = nd, type = "prob",
                         level = "subject"), regexp = "[Uu]nknown|population")
  })
  it("glance reports mode=descriptive and n_random_effects=2", {
    g <- glance(.get_desc_fit())
    expect_equal(g$mode, "descriptive")
    expect_equal(g$n_random_effects, 2L)
  })
  it("augment + summary + print + print(summary) run for descriptive", {
    au <- augment(.get_desc_fit())
    expect_true(all(c(".fitted", ".resid") %in% names(au)))
    expect_s3_class(summary(.get_desc_fit()), "summary.beezdiscounting_choice")
    expect_invisible(print(.get_desc_fit()))
    expect_invisible(print(summary(.get_desc_fit())))
  })
  it("confint has distinct RE-SD labels and a bounded (tanh) correlation CI", {
    ci <- suppressWarnings(confint(.get_desc_fit()))
    expect_false(any(duplicated(ci$term)))             # no duplicate row labels
    expect_true("cor_slopes" %in% ci$term)
    cr <- ci[ci$term == "cor_slopes", ]
    expect_true(cr$estimate >= -1 && cr$estimate <= 1) # correlation-scaled (tanh)
    # CI bounds lie in [-1, 1] when SEs are reliable; NA otherwise (degenerate fit)
    if (is.finite(cr$conf.low) && is.finite(cr$conf.high)) {
      expect_true(cr$conf.low >= -1 && cr$conf.high <= 1)
    }
  })
})

describe("VarCorr (structural) still reports sigma_u", {
  it("structural VarCorr surfaces the log-k RE SD", {
    dat <- .choice_fit_fixture(n_subjects = 30, seed = 41)
    fit <- fit_dd_choice(dat, mode = "structural", verbose = 0)
    vc <- nlme::VarCorr(fit)
    expect_true(is.data.frame(vc) || is.matrix(vc) || is.list(vc))
  })
})

describe("descriptive S3 — pooled (random_slopes = FALSE, q = 0)", {
  .pooled <- new.env(parent = emptyenv())
  get_pooled <- function() {
    if (!exists("fit", envir = .pooled)) {
      dat <- .choice_desc_fixture(n_subjects = 50, seed = 33)
      .pooled$fit <- fit_dd_choice(dat, mode = "descriptive",
                                   random_slopes = FALSE, verbose = 0)
    }
    .pooled$fit
  }
  it("ranef returns an id-only frame (no slopes)", {
    re <- ranef(get_pooled())
    expect_true("id" %in% names(re))
    expect_false(any(c("b_mag", "b_delay") %in% names(re)))
  })
  it("VarCorr returns an empty (no-RE) structure", {
    vc <- nlme::VarCorr(get_pooled())
    expect_equal(nrow(vc), 0L)
  })
  it("predict(type='prob') works; subject == population (no RE)", {
    p_sub <- predict(get_pooled(), type = "prob", level = "subject")$.prob
    p_pop <- predict(get_pooled(), type = "prob", level = "population")$.prob
    expect_equal(p_sub, p_pop, tolerance = 1e-10)
    expect_true(all(p_sub > 0 & p_sub < 1))
  })
  it("predict(type='parameters') returns the id-only slope table", {
    pp <- predict(get_pooled(), type = "parameters")
    expect_true("id" %in% names(pp))
  })
  it("glance reports n_random_effects = 0", {
    expect_equal(glance(get_pooled())$n_random_effects, 0L)
  })
})

describe("descriptive RE extraction orientation", {
  it("ranef slope rows are in subject_levels order", {
    fit <- .get_desc_fit()
    re <- ranef(fit)
    expect_identical(as.character(re$id),
                     as.character(fit$param_info$subject_levels))
  })
})
