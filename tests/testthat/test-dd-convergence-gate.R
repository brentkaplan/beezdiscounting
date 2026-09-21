# Audit 2026-09-06 F-BZ4-1: multi-start selection must prefer converged starts,
# non-convergence must be announced at fit time regardless of `verbose`, and the
# Wald surfaces (tidy / confint / summary / EMMs / comparisons) must warn on a
# non-converged fit instead of presenting its p-values and CIs silently.

cand <- function(nll, conv, idx) {
  list(nll = nll, opt = list(convergence = conv), start_idx = idx)
}
no_blowup <- function(opt) FALSE

describe(".dd_select_start()", {
  it("prefers a converged start over a lower-NLL non-converged one", {
    best <- .dd_select_start(
      list(cand(10, 1L, 1L), cand(12, 0L, 2L), cand(15, 0L, 3L)), no_blowup)
    expect_identical(best$start_idx, 2L)
    info <- best$selection
    expect_identical(info$tier, 1L)
    expect_identical(info$n_converged, 2L)
    expect_true(info$lower_nll_nonconverged)
  })

  it("falls back to the lowest-NLL non-blowup start when none converged", {
    best <- .dd_select_start(
      list(cand(10, 1L, 1L), cand(8, 1L, 2L), NULL), no_blowup)
    expect_identical(best$start_idx, 2L)
    expect_identical(best$selection$tier, 2L)
    expect_identical(best$selection$n_starts, 3L)
    expect_identical(best$selection$n_finite, 2L)
  })

  it("does not accept a converged blow-up over a non-converged sane start", {
    blow <- function(opt) isTRUE(opt$blow)
    a <- cand(5, 0L, 1L); a$opt$blow <- TRUE
    b <- cand(9, 1L, 2L)
    best <- .dd_select_start(list(a, b), blow)
    expect_identical(best$start_idx, 2L)
    expect_identical(best$selection$tier, 2L)
  })

  it("returns tier 3 when every start blew up and NULL when none is finite", {
    blow <- function(opt) TRUE
    best <- .dd_select_start(list(cand(5, 0L, 1L), cand(3, 1L, 2L)), blow)
    expect_identical(best$start_idx, 2L)
    expect_identical(best$selection$tier, 3L)
    expect_null(.dd_select_start(list(NULL, cand(Inf, 0L, 2L)), no_blowup))
  })

  it("reports no lower-NLL non-converged start when the winner is the minimum", {
    best <- .dd_select_start(list(cand(10, 0L, 1L), cand(12, 1L, 2L)), no_blowup)
    expect_false(best$selection$lower_nll_nonconverged)
  })
})

describe("fit_dd_tmb() non-convergence gate", {
  sim <- simulate_dd_ip(n_subjects = 30, seed = 3)

  it("warns at fit time with verbose = 0 and records the selection", {
    expect_warning(
      fit <- fit_dd_tmb(sim, equation = "mazur", verbose = 0,
                        tmb_control = list(iter_max = 3)),
      class = "beezdiscounting_convergence_warning"
    )
    expect_false(fit$converged)
    expect_identical(fit$multi_start_info$n_converged, 0L)
    expect_identical(fit$multi_start_info$tier, 2L)
    expect_warning(tidy(fit), class = "beezdiscounting_convergence_warning")
    expect_warning(confint(fit), class = "beezdiscounting_convergence_warning")
    expect_warning(summary(fit), class = "beezdiscounting_convergence_warning")
  })

  it("warns on the single-start path too", {
    expect_warning(
      fit_dd_tmb(sim, equation = "mazur", verbose = 0, multi_start = FALSE,
                 tmb_control = list(iter_max = 3)),
      class = "beezdiscounting_convergence_warning"
    )
  })

  it("stays silent for a converged fit", {
    fit <- expect_no_warning(fit_dd_tmb(sim, equation = "mazur", verbose = 0))
    expect_true(fit$converged)
    expect_identical(fit$multi_start_info$tier, 1L)
    expect_no_warning(tidy(fit))
    expect_no_warning(confint(fit))
  })

  it("warns from the EMM and comparison surfaces", {
    sim2 <- simulate_dd_ip(n_subjects = 30, n_conditions = 2,
                           delta_k = c(0, log(2)), seed = 3)
    fit <- suppressWarnings(
      fit_dd_tmb(sim2, equation = "mazur", factors = "condition", verbose = 0,
                 tmb_control = list(iter_max = 3)))
    expect_false(fit$converged)
    expect_warning(get_dd_param_emms(fit, factors_in_emm = "condition"),
                   class = "beezdiscounting_convergence_warning")
    expect_warning(get_dd_comparisons(fit),
                   class = "beezdiscounting_convergence_warning")
  })
})

describe("fit_dd_choice() non-convergence gate", {
  it("structural: warns at fit time and from tidy/confint", {
    dat <- .choice_fit_fixture(n_subjects = 30)
    expect_warning(
      fit <- fit_dd_choice(dat, mode = "structural", verbose = 0,
                           tmb_control = list(iter_max = 2)),
      class = "beezdiscounting_convergence_warning"
    )
    expect_false(fit$converged)
    expect_false(is.null(fit$multi_start_info))
    expect_warning(tidy(fit), class = "beezdiscounting_convergence_warning")
    expect_warning(confint(fit), class = "beezdiscounting_convergence_warning")
  })

  it("descriptive: warns at fit time and from tidy/confint", {
    dat <- .choice_desc_fixture(n_subjects = 30)
    # the truncated fit also trips sdreport (NaN / non-PD); collect every class
    classes <- character(0)
    fit <- withCallingHandlers(
      fit_dd_choice(dat, mode = "descriptive", verbose = 0,
                    tmb_control = list(iter_max = 2)),
      warning = function(w) {
        classes <<- c(classes, class(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_true("beezdiscounting_convergence_warning" %in% classes)
    expect_false(fit$converged)
    expect_warning(tidy(fit), class = "beezdiscounting_convergence_warning")
    expect_warning(
      expect_warning(confint(fit), class = "beezdiscounting_convergence_warning"),
      "unreliable"
    )
  })
})
