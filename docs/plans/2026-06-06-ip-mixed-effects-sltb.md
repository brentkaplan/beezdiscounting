# beezdiscounting IP Mixed-Effects (SLT-beta + Gaussian, TMB) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a TMB-based mixed-effects discounting tier to beezdiscounting — `fit_dd_tmb()` fitting indifference-point data with a subject random intercept on `log k` under a Scale-Location-Truncated Beta (default) or Gaussian likelihood, plus the full S3 + emmeans surface, simulation, and a user vignette.

**Architecture:** Port beezdemand's TMB tier as the single-random-intercept special case. The effects sit on `log k = Xβ + σ·u` (linear in β), so EMMs/contrasts are computed for `k` on the log scale and back-transformed (the "EMMs-without-emmeans" averaging-matrix trick) — *family-agnostic*. The observation `family` (`sltb`/`gaussian`) and `equation` (`mazur`/`exponential`) are integer dispatch dimensions in one C++ template (`src/MixedDiscounting.h`); the family only swaps the density term and the auxiliary parameter (`log_aux` = `log φ` or `log σ_e`).

**Tech Stack:** R, TMB (>= 1.9.0) + RcppEigen (compiled C++), emmeans/broom/generics contracts, testthat 3e (BDD describe/it), roxygen2, knitr vignette, GitHub Actions CI.

**Design spec:** `docs/specs/2026-06-06-ip-mixed-effects-sltb-design.md` (read it first). **Verified reference density:** `dev/sltb-verification/` (the SLT density is confirmed correct; constants `s=1.0000001`, `l=1e-8`; `−log(Z)` load-bearing).

---

## Integration & execution notes (READ FIRST — resolves cross-phase contracts)

The plan is grouped into phases F (Foundation, pure R) → T (TMB core) → P (Fit pipeline) → M (S3 methods) → E (emmeans) → S (Surface). Tasks are id'd `Task <PHASE>.<n>`. The following reconciliations override any drift in individual phase sections:

1. **Density helpers live in `R/dd-density.R`**: `.dd_slt_logpdf(y, mu, phi, s=1.0000001, l=1e-8)` and `.dd_gaussian_logpdf(y, mu, sigma_e)` (Task F.1). The compile-gate test (T.5) and any method needing the R density **source these** — ignore T.5's inline `.gate_slt_logpdf` fallback (keep only if you want belt-and-suspenders; prefer the real symbol).

2. **Single RHS/design builder (load-bearing for emmeans correctness).** Port `build_fixed_rhs(factors, factor_interaction, continuous_covariates, data)` from beezdemand `R/utils.R` into beezdiscounting `R/utils.R` as a shared helper (add this as **Task P.2a**, executed inside Phase P before P.2 finishes). `.dd_tmb_build_design` (P.2) uses it and **stores both the `rhs` formula and `attr(X, "contrasts")`** in `fit$formula_details` (`formula_details = list(X, rhs, contrasts)`). `.dd_build_emm_ref_grid` (E.1) rebuilds the grid design with `model.matrix(rhs, grid, contrasts.arg = fit$formula_details$contrasts)` and reorders to the fitted column order, aborting on mismatch. This guarantees EMM columns align with the fitted `beta_k` — do not let P and E construct designs by different routes.

3. **Coefficient naming.** `fit$model$coefficients` keeps the optimizer names: the fixed-effect block is named `beta_k` (length `ncol(X)`; when `ncol(X) > 1`, R names them `beta_k1, beta_k2, …`). The **population log-k intercept** is the entry for X's intercept column (the first `beta_k`). Reference the block via `grepl("^beta_k", names(...))` and the intercept via its first element — never assume a scalar `coefficients[["beta_k"]]` when factors are present. M's display term names (`k:(Intercept)`, `k:groupB`) are derived in `tidy`/`confint` only and must not rename the stored coefficients. Tests in P.8 and S.2 use the first `beta_k` element as the population intercept.

4. **`log_aux` rename.** The optimizer parameter is `log_aux`; `.dd_tmb_extract_estimates` (P.6) renames it to `log_phi` (sltb) / `log_sigma_e` (gaussian) in `model$coefficients` and aligned `se`. All downstream (M, E, S) use the renamed form; `opt$par` retains `log_aux`.

5. **Accepted signature deltas from the contract (intentional):** `.dd_tmb_default_starts(prepared, design, family, equation="mazur")` (extra `equation` for data-driven start inversion); `.dd_tmb_run_optimizer(obj, start, tmb_control, user_specified, verbose)` and `.dd_tmb_multi_start(tmb_data, starts_list, tmb_control, user_specified, verbose)` (extra `user_specified` from the verbatim beezdemand port). These are fine; keep them consistent across P.

6. **One `%||%`.** Define `%||%` once (in `R/utils.R`, or reuse rlang's). Phases that wrote a local `%||%` must drop it in favor of the single package definition (dedupe at integration).

7. **DESCRIPTION edited once (Task T.1).** T.1 is the canonical DESCRIPTION task: add `LinkingTo: TMB, RcppEigen`; `Imports: TMB (>= 1.9.0), RcppEigen, emmeans, generics`; `Suggests: knitr, rmarkdown, testthat`; `VignetteBuilder: knitr`; `SystemRequirements: GNU make`; `Version: 0.4.0`; `Date: 2026-06-06`. S.7 must **not** re-edit DESCRIPTION (only NEWS/`devtools::document()`/CI). Note `generics` (needed by E.4 + M for tidy/glance/augment re-exports) is added here.

8. **One NEWS.md.** F.3 creates `NEWS.md` with the `# beezdiscounting 0.4.0` header; every later phase **prepends bullets under that same header** — never recreate the file. S.7 finalizes the 0.4.0 entry.

9. **Execution order:** F.1–F.4 → T.1–T.5 → P.2a (build_fixed_rhs) → P.1–P.8 → **S.1 (simulator + `helper-dd-sim.R`)** → M.0–M.8 → E.1–E.4 → S.2–S.7. The simulator (S.1) is pulled forward because M and E tests depend on `.simulate_dd_ip_mixed`. Commit after every task. Each compiled-template test starts with `skip_on_cran()` + `skip_if_not_installed("TMB")`; pure-R tests need no skips.

10. **Subject-par `phi`:** MVP keeps `φ` population-level, so `ranef`/`subject_pars` have no `phi` column (M asserts its absence). The degenerate φ→0/k→∞ optimum (spec §4.8) is guarded in `.dd_tmb_multi_start` (P.5, `.dd_phi_min = 0.1`) with a regression test in P.8.

---
## Phase F — Foundation (pure R, no compile)

This phase builds the pure-R substrate that the compiled TMB template (Phase G) is later
cross-checked against, plus the parameter-space machinery, the input validator, and the
random-effects normalizer. Nothing here requires a C++ toolchain, so it is the fastest first
slice (spec §7.1–§7.2). All work happens on branch `feat/tmb-mixed-discounting` in
`/Users/brent/Dropbox/GIT/beezdiscounting`. Tests use testthat 3e with BDD `describe()`/`it()`.

### Task F.1: SLT-beta + Gaussian log-density R helpers (the source of truth)

The `.dd_slt_logpdf` function is transcribed verbatim from the VERIFIED reference
`dev/sltb-verification/verify_sltb.R::slt_logpdf` (constants `s = 1.0000001`, `l = 1e-8`).
This R function is the canonical density that the C++ template's `-nll` is later cross-checked
against to 1e-8 (Phase G). `.dd_gaussian_logpdf` is the baseline comparator family.

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-density.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-sltb-density.R`

- [ ] **Step 1: Write failing test `test-sltb-density.R`** (formalizes the six property
  checks from `dev/sltb-verification/verify_sltb.R`: kernel↔Z normalization, normalized
  density integrates to 1, Z load-bearing, boundary finiteness at 0/1, beta-limit `f=dbeta/Z`,
  exact match to the reference body, plus the Gaussian helper).

```r
# tests/testthat/test-sltb-density.R
# Property tests for the pure-R SLT-beta + Gaussian log-densities.
# These formalize dev/sltb-verification/verify_sltb.R as package tests; the
# C++ template (-nll) is later cross-checked against .dd_slt_logpdf to 1e-8.

S_REF <- 1.0000001
L_REF <- 1e-8

# Analytic truncation normalizer Z (matches verify_sltb.R::Zof)
Zof <- function(mu, phi, s = S_REF, l = L_REF) {
  a <- mu * phi
  b <- (1 - mu) * phi
  pbeta(1 / s + l, a, b) - pbeta(l, a, b)
}

# Exact reference function body from verify_sltb.R (sum of -ll over rows)
ref_nll_mazur <- function(data, par) {
  k <- par[1]
  phi <- par[2]
  delay <- data$delay
  ip <- data$IP
  mu <- 1 / (k * delay + 1)
  alpha <- mu * phi
  beta <- phi * (1 - mu)
  s <- 1.0000001
  l <- 0.00000001
  ll_temp <- lgamma(alpha + beta) - lgamma(beta) - lgamma(alpha) +
    (alpha - 1) * log(ip / s + l) + (beta - 1) * log(1 - (ip / s + l)) -
    log(s) - log(pbeta(1 / s + l, alpha, beta) - pbeta(l, alpha, beta))
  sum(-(ll_temp * ifelse(ip >= 0 & ip <= 1, 1, 0)))
}

describe(".dd_slt_logpdf", {
  mu_grid <- c(0.02, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99)
  phi_grid <- c(1.5, 2, 5, 20, 100, 500)
  grid <- expand.grid(mu = mu_grid, phi = phi_grid)

  it("normalizes: kernel integral equals analytic Z (moderate shapes)", {
    mod <- with(grid, mu * phi >= 1 & (1 - mu) * phi >= 1)
    kern_err <- mapply(function(mu, phi) {
      a <- mu * phi
      b <- (1 - mu) * phi
      num <- stats::integrate(
        function(g) (1 / S_REF) * dbeta(g / S_REF + L_REF, a, b),
        0, 1, rel.tol = 1e-10, subdivisions = 1000L
      )$value
      abs(num - Zof(mu, phi))
    }, grid$mu[mod], grid$phi[mod])
    expect_lt(max(kern_err), 1e-7)
  })

  it("integrates the normalized density to 1 (moderate shapes)", {
    mod <- with(grid, mu * phi >= 1 & (1 - mu) * phi >= 1)
    norm_err <- mapply(function(mu, phi) {
      abs(stats::integrate(
        function(g) exp(.dd_slt_logpdf(g, mu, phi)),
        0, 1, rel.tol = 1e-10
      )$value - 1)
    }, grid$mu[mod], grid$phi[mod])
    expect_lt(max(norm_err), 1e-6)
  })

  it("keeps Z load-bearing (Z varies with shape, not always 1)", {
    expect_lt(Zof(0.05, 1), 0.7)
    expect_gt(Zof(0.5, 10), 0.99)
  })

  it("is finite at the boundaries y=0 and y=1 for all shapes", {
    b0 <- mapply(function(mu, phi) .dd_slt_logpdf(0, mu, phi), grid$mu, grid$phi)
    b1 <- mapply(function(mu, phi) .dd_slt_logpdf(1, mu, phi), grid$mu, grid$phi)
    expect_true(all(is.finite(c(b0, b1))))
  })

  it("matches dbeta/Z in the interior (scale/location limit)", {
    beta_lim_err <- mapply(function(mu, phi) {
      ips <- c(0.1, 0.3, 0.5, 0.7, 0.9)
      max(abs(exp(.dd_slt_logpdf(ips, mu, phi)) -
        dbeta(ips, mu * phi, (1 - mu) * phi) / Zof(mu, phi)))
    }, grid$mu, grid$phi)
    expect_lt(max(beta_lim_err), 1e-3)
  })

  it("equals the reference NLL body exactly (incl. boundary IP=0)", {
    dtest <- data.frame(
      delay = c(7, 30, 180, 365, 730, 1460, 2920),
      IP = c(0.95, 0.8, 0.5, 0.35, 0.2, 0.08, 0.0)
    )
    k0 <- 0.01
    phi0 <- 8
    mu_m <- 1 / (1 + k0 * dtest$delay)
    ours <- sum(-.dd_slt_logpdf(dtest$IP, mu_m, phi0))
    expect_equal(ours, ref_nll_mazur(dtest, c(k0, phi0)), tolerance = 1e-12)
  })

  it("is vectorized over y, mu, and phi", {
    y <- c(0, 0.5, 1)
    mu <- c(0.2, 0.5, 0.8)
    phi <- c(5, 10, 20)
    out <- .dd_slt_logpdf(y, mu, phi)
    expect_length(out, 3L)
    expect_true(all(is.finite(out)))
  })
})

describe(".dd_gaussian_logpdf", {
  it("equals dnorm(log=TRUE)", {
    y <- c(0, 0.25, 0.5, 0.75, 1)
    mu <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    sigma_e <- 0.1
    expect_equal(
      .dd_gaussian_logpdf(y, mu, sigma_e),
      dnorm(y, mu, sigma_e, log = TRUE),
      tolerance = 1e-12
    )
  })

  it("is vectorized and finite for positive sigma_e", {
    out <- .dd_gaussian_logpdf(c(0.2, 0.8), c(0.3, 0.7), 0.05)
    expect_length(out, 2L)
    expect_true(all(is.finite(out)))
  })
})
```

- [ ] **Step 2: Run the test, expect FAIL** (functions do not exist yet).

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-sltb-density.R")'
```
Expected output: errors like `could not find function ".dd_slt_logpdf"` / `".dd_gaussian_logpdf"`; test run reports failures (non-zero).

- [ ] **Step 3: Implement `R/dd-density.R`** (verbatim port of the verified reference; note
  the lgamma sign order `lgamma(a+b) - lgamma(b) - lgamma(a)` is algebraically identical to
  the spec's `lgamma(a+b) - lgamma(a) - lgamma(b)` — addition commutes — and is copied exactly
  from the reference for byte-for-byte tie-out).

```r
# R/dd-density.R

#' SLT-beta and Gaussian log-densities (pure R reference)
#'
#' Pure-R implementations of the two observation-family log-densities used by
#' the mixed-effects discounting model. `.dd_slt_logpdf()` is transcribed
#' verbatim from the verified reference
#' `dev/sltb-verification/verify_sltb.R::slt_logpdf` and is the source of truth
#' against which the compiled TMB template (`src/MixedDiscounting.h`) is
#' cross-checked to 1e-8. See Kim, Koffarnus & Franck (2024) and
#' Kim, Kaplan, Koffarnus & Franck (2025; arXiv:2509.13167).
#'
#' The scale-location-truncation (SLT) constants are the *reference-code*
#' values `s = 1.0000001` (= 1 + 1e-7) and `l = 1e-8`. The truncation
#' normalizer `Z = pbeta(1/s + l, a, b) - pbeta(l, a, b)` is load-bearing
#' (NOT approx 1 at small shapes) and the boundaries are finite at y = 0 and
#' y = 1 by construction.
#'
#' @param y Numeric vector of indifference proportions in `[0, 1]`.
#' @param mu Numeric (scalar or vector recycled against `y`); the mean on the
#'   identity link (a discounting function value).
#' @param phi Numeric (scalar or vector); SLT-beta precision (> 0).
#' @param s,l SLT scale/location constants; defaults are the verified values.
#' @return Numeric vector of log-density values, the same length as the
#'   recycled inputs.
#' @keywords internal
#' @noRd
.dd_slt_logpdf <- function(y, mu, phi, s = 1.0000001, l = 1e-8) {
  a <- mu * phi
  b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(b) - lgamma(a) +
    (a - 1) * log(y / s + l) + (b - 1) * log(1 - (y / s + l)) -
    log(s) - log(pbeta(1 / s + l, a, b) - pbeta(l, a, b))
}

#' @rdname dot-dd_slt_logpdf
#' @param sigma_e Numeric (> 0); Gaussian residual standard deviation.
#' @keywords internal
#' @noRd
.dd_gaussian_logpdf <- function(y, mu, sigma_e) {
  stats::dnorm(y, mean = mu, sd = sigma_e, log = TRUE)
}
```

- [ ] **Step 4: Run the test, expect PASS.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-sltb-density.R")'
```
Expected output: all `it()` blocks pass; `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 9 ]`.

- [ ] **Step 5: Commit (targeted add).**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add R/dd-density.R tests/testthat/test-sltb-density.R && \
git commit -m "feat(density): add pure-R SLT-beta and Gaussian log-densities

Transcribe .dd_slt_logpdf verbatim from the verified
dev/sltb-verification/verify_sltb.R (s=1.0000001, l=1e-8); this is the
source of truth the C++ template is later cross-checked against to 1e-8.
Add .dd_gaussian_logpdf baseline and formalize the six density property
checks (normalization, Z, boundary finiteness, beta-limit, ref-body match)
as test-sltb-density.R.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task F.2: Parameter registry + param-space delta-method transforms

Port `.dd_param_registry`, `.dd_transform_est_se`, and `.dd_transform_coef_table` from
beezdemand. The registry keys are `"k"`, `"phi"`, `"s"` (per the pinned contract). The
core-term regex is retargeted from beezdemand's demand terms to
`^k($|_)|^s($|_)|^phi($|_)`.

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-param-registry.R`
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-param-space.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-param-space.R`

- [ ] **Step 1: Write failing test `test-dd-param-space.R`** (delta-method round-trips,
  registry shape, term-display regex, coef-table retargeting).

```r
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
    # 'kappa_not_core' must NOT match ^k($|_) (next char is 'a', not $ or _)
    expect_equal(out$estimate[out$term == "kappa_not_core"], 5)
  })
})
```

- [ ] **Step 2: Run the test, expect FAIL** (registry + functions absent).

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-dd-param-space.R")'
```
Expected output: `object '.dd_param_registry' not found` / `could not find function ".dd_transform_est_se"`; failures reported.

- [ ] **Step 3: Implement `R/dd-param-registry.R`** (new registry keyed by k/phi/s; mirrors
  the field shape of `beezdemand/R/param-registry.R:78-105` `.beezdemand_param_registry`,
  reduced to the three discounting parameters from the pinned contract).

```r
# R/dd-param-registry.R

#' Parameter naming registry for beezdiscounting mixed-effects models
#'
#' Single source of truth for the canonical parameters of the IP-family
#' mixed-effects discounting model: the subject discount rate `k`, the
#' SLT-beta precision `phi`, and the (currently fixed) scale constant `s`.
#' Fields mirror `beezdemand`'s `.beezdemand_param_registry` (canonical,
#' description, constraint, valid_scales, default_scale).
#'
#' @keywords internal
#' @noRd
.dd_param_registry <- list(
  k = list(
    canonical = "k",
    description = "Discount rate (Mazur hyperbolic / exponential)",
    constraint = "k > 0",
    valid_scales = c("natural", "log", "log10"),
    default_scale = "natural"
  ),

  phi = list(
    canonical = "phi",
    description = "SLT-beta precision (population)",
    constraint = "phi > 0",
    valid_scales = c("natural", "log", "log10"),
    default_scale = "natural"
  ),

  s = list(
    canonical = "s",
    description = "SLT scale constant (fixed at 1.0000001 in the MVP)",
    constraint = "s >= 1",
    valid_scales = c("natural", "log", "log10"),
    default_scale = "natural"
  )
)
```

- [ ] **Step 4: Implement `R/dd-param-space.R`** (port `beezdemand_transform_est_se`
  =`beezdemand/R/param-space.R:115-178` and `beezdemand_transform_coef_table`
  =`beezdemand/R/param-space.R:180-237`, plus the supporting validators
  `:20-37` and the term-display helper. Adaptations: rename `beezdemand_` → `.dd_`; the core
  regex in `.dd_transform_coef_table` becomes `^k($|_)|^s($|_)|^phi($|_)`; `.dd_term_display_space`
  handles only k/s/phi prefixes; the `%||%` null-coalescer is defined locally to avoid a
  cross-package import).

```r
# R/dd-param-space.R

# Local null-coalescing operator (ports the rlang %||% used by beezdemand's
# param-space helpers) so this file has no extra import.
`%|0|%` <- function(x, y) if (is.null(x)) y else x

#' Validate a report-space string
#'
#' Ported from beezdemand::beezdemand_validate_report_space
#' (R/param-space.R:20-37).
#' @keywords internal
#' @noRd
.dd_validate_report_space <- function(
  report_space,
  choices = c("natural", "log", "log10", "internal")
) {
  if (is.null(report_space) || length(report_space) != 1) {
    stop("'report_space' must be a single character value.", call. = FALSE)
  }
  report_space <- as.character(report_space)
  if (!report_space %in% choices) {
    stop(
      "'report_space' must be one of: ",
      paste(sprintf('"%s"', choices), collapse = ", "), ".",
      call. = FALSE
    )
  }
  report_space
}

#' Display label for a core term in a given report space
#'
#' Retargeted from beezdemand::beezdemand_term_display_space
#' (R/param-space.R:39-90) to the discounting core terms k / s / phi.
#' @keywords internal
#' @noRd
.dd_term_display_space <- function(term, report_space) {
  if (is.na(term) || is.null(term)) return(NA_character_)
  term <- as.character(term)
  report_space <- as.character(report_space)

  label <- function(base, suffix) {
    prefix <- switch(report_space,
      log10 = sprintf("log10(%s)", base),
      log = sprintf("log(%s)", base),
      base
    )
    paste0(prefix, suffix)
  }

  if (grepl("^k($|_)", term)) return(label("k", sub("^k", "", term)))
  if (grepl("^phi($|_)", term)) return(label("phi", sub("^phi", "", term)))
  if (grepl("^s($|_)", term)) return(label("s", sub("^s", "", term)))
  term
}

#' Delta-method transform of an estimate + SE between parameter spaces
#'
#' Verbatim port of beezdemand::beezdemand_transform_est_se
#' (R/param-space.R:115-178); only the function name changes.
#'
#' @param estimate,se Numeric estimate and its standard error.
#' @param from,to One of "natural", "log", "log10", "internal".
#' @return list(estimate=, se=).
#' @keywords internal
#' @noRd
.dd_transform_est_se <- function(estimate, se, from, to) {
  if (identical(from, to) || to == "internal") {
    return(list(estimate = estimate, se = se))
  }

  ln10 <- log(10)

  if (from == "natural" && to == "log10") {
    new_est <- ifelse(is.finite(estimate) & estimate > 0, log10(estimate), NA_real_)
    new_se <- ifelse(is.finite(se) & is.finite(estimate) & estimate > 0,
      se / (estimate * ln10), NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log10" && to == "natural") {
    new_est <- 10^estimate
    new_se <- ifelse(is.finite(se), ln10 * (10^estimate) * se, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log" && to == "natural") {
    new_est <- exp(estimate)
    new_se <- ifelse(is.finite(se), exp(estimate) * se, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "natural" && to == "log") {
    new_est <- ifelse(is.finite(estimate) & estimate > 0, log(estimate), NA_real_)
    new_se <- ifelse(is.finite(se) & is.finite(estimate) & estimate > 0,
      se / estimate, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log" && to == "log10") {
    new_est <- estimate / ln10
    new_se <- ifelse(is.finite(se), se / ln10, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log10" && to == "log") {
    new_est <- estimate * ln10
    new_se <- ifelse(is.finite(se), se * ln10, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }

  stop("Unsupported transform from '", from, "' to '", to, "'.", call. = FALSE)
}

#' Transform the core (k/phi/s) rows of a coefficient table to a report space
#'
#' Port of beezdemand::beezdemand_transform_coef_table
#' (R/param-space.R:180-237). The core-term detector regex is retargeted to
#' `^k($|_)|^s($|_)|^phi($|_)`; everything else (column bookkeeping, the
#' per-row delta-method call, display labels) is unchanged.
#'
#' @keywords internal
#' @noRd
.dd_transform_coef_table <- function(
  coef_tbl,
  report_space,
  internal_space,
  term_col = "term",
  estimate_col = "estimate",
  se_col = "std.error",
  scale_col = "estimate_scale",
  display_col = "term_display"
) {
  report_space <- .dd_validate_report_space(report_space)

  if (!nrow(coef_tbl)) return(coef_tbl)
  if (!all(c(term_col, estimate_col, se_col) %in% names(coef_tbl))) {
    return(coef_tbl)
  }

  is_core <- grepl("^k($|_)", coef_tbl[[term_col]]) |
    grepl("^s($|_)", coef_tbl[[term_col]]) |
    grepl("^phi($|_)", coef_tbl[[term_col]])

  out <- coef_tbl
  if (!("estimate_internal" %in% names(out))) out$estimate_internal <- NA_real_
  if (!(display_col %in% names(out))) out[[display_col]] <- as.character(out[[term_col]])
  if (!(scale_col %in% names(out))) out[[scale_col]] <- NA_character_

  for (i in which(is_core)) {
    term <- as.character(out[[term_col]][i])

    from_space <- out[[scale_col]][i] %|0|% internal_space
    if (is.na(from_space) || !nzchar(from_space)) from_space <- internal_space

    to_space <- report_space
    if (to_space == "internal") to_space <- from_space

    trans <- .dd_transform_est_se(
      estimate = out[[estimate_col]][i],
      se = out[[se_col]][i],
      from = from_space,
      to = to_space
    )

    out$estimate_internal[i] <- out[[estimate_col]][i]
    out[[estimate_col]][i] <- trans$estimate
    out[[se_col]][i] <- trans$se

    out[[scale_col]][i] <- to_space
    out[[display_col]][i] <- .dd_term_display_space(term, to_space)
  }

  out
}
```

- [ ] **Step 5: Run the test, expect PASS.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-dd-param-space.R")'
```
Expected output: all blocks pass; `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 12 ]`.

- [ ] **Step 6: Commit (targeted add).**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add R/dd-param-registry.R R/dd-param-space.R tests/testthat/test-dd-param-space.R && \
git commit -m "feat(param-space): port registry + delta-method transforms for k/phi/s

Add .dd_param_registry (keyed by k, phi, s) and port
.dd_transform_est_se / .dd_transform_coef_table from
beezdemand/R/param-space.R, retargeting the core-term regex to
^k(\$|_)|^s(\$|_)|^phi(\$|_). Covers natural<->log<->log10 delta-method
round-trips and core-row selection.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task F.3: `.dd_validate_ip` — coerce / clamp / WARN loudly + document

The IP validator detects percent/amount-scaled inputs and divides them to `[0,1]`, clamps
mild out-of-range values, and WARNS loudly (naming counts), per spec §5. It returns the long
`id, x, y` data frame plus a `coercion_info` audit list. It is the single coercion choke-point
so the likelihood never sees out-of-range `y` (spec §4.1: filter, never mask).

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-validate.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-validate.R`

- [ ] **Step 1: Write failing test `test-dd-validate.R`** (proportion pass-through, percent
  ÷100 with warning, amount ÷ll with warning, overshoot clamp with named-count warning, NA-y
  error, missing-column error, column remap).

```r
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
      regexp = "clamp"
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
})
```

- [ ] **Step 2: Run the test, expect FAIL** (validator absent).

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-dd-validate.R")'
```
Expected output: `could not find function ".dd_validate_ip"`; failures reported.

- [ ] **Step 3: Implement `R/dd-validate.R`** (complete; loud `warning()` with named counts,
  documented behavior in roxygen per spec §5).

```r
# R/dd-validate.R

#' Validate and coerce IP-family long data to proportions in [0, 1]
#'
#' The single coercion choke-point for indifference-point (IP) mixed-effects
#' discounting. It remaps caller column names to the canonical `id, x, y` long
#' format (matching `fit_dd()`), detects and divides percent- or amount-scaled
#' responses to the `[0, 1]` proportion scale, and clamps mild post-coercion
#' overshoot. **All coercion is loud**: percent/amount division and any
#' clamping emit a `warning()` that names how many values were affected. This
#' guarantees the downstream likelihood never sees out-of-range `y` (we filter
#' and coerce here rather than masking inside the likelihood; see the design
#' spec, "Implementation landmines").
#'
#' Scale detection:
#' \itemize{
#'   \item `response_scale = "proportion"` (default): if `max(y, na.rm) > 1.5`
#'     the data are treated as percent and divided by 100 (with a warning);
#'     otherwise passed through.
#'   \item `response_scale = "percent"`: always divide by 100.
#'   \item `response_scale = "amount"`: divide by `ll` (the larger-later
#'     reward); `ll` is required.
#' }
#' After scaling, values `> 1` are clamped to 1 and values `< 0` to 0, each
#' with a warning naming the count. `0` and `1` are valid (the SLT-beta
#' family's purpose), so they are never warned about.
#'
#' @param data A data frame containing the id, delay (x), and indifference
#'   proportion (y) columns.
#' @param y_var,x_var,id_var Column names in `data` for the response, delay,
#'   and subject id. Defaults `"y"`, `"x"`, `"id"`.
#' @param ll Optional numeric larger-later reward; required when
#'   `response_scale = "amount"`.
#' @param response_scale One of `"proportion"` (default), `"percent"`,
#'   `"amount"`.
#' @return A list with:
#'   \describe{
#'     \item{data}{a data frame with exactly columns `id`, `x`, `y`.}
#'     \item{coercion_info}{list(divided_by, n_clamped_hi, n_clamped_lo,
#'       scale_detected).}
#'   }
#' @keywords internal
#' @noRd
.dd_validate_ip <- function(data,
                            y_var = "y",
                            x_var = "x",
                            id_var = "id",
                            ll = NULL,
                            response_scale = c("proportion", "percent", "amount")) {
  response_scale <- match.arg(response_scale)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  needed <- c(id_var, x_var, y_var)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "Column(s) not found in `data`: ",
      paste(shQuote(missing_cols), collapse = ", "),
      ". Set id_var / x_var / y_var to match your data.",
      call. = FALSE
    )
  }

  long <- data.frame(
    id = data[[id_var]],
    x = as.numeric(data[[x_var]]),
    y = as.numeric(data[[y_var]]),
    stringsAsFactors = FALSE
  )

  # --- scale detection / division ------------------------------------------
  max_y <- suppressWarnings(max(long$y, na.rm = TRUE))
  divided_by <- 1
  scale_detected <- "proportion"

  if (response_scale == "amount") {
    if (is.null(ll) || !is.finite(ll) || ll <= 0) {
      stop(
        "response_scale = 'amount' requires a positive larger-later reward `ll`.",
        call. = FALSE
      )
    }
    divided_by <- ll
    scale_detected <- "amount"
    long$y <- long$y / ll
    warning(
      sprintf(
        "Detected amount-scale responses; divided y by the larger-later reward (ll = %s).",
        format(ll)
      ),
      call. = FALSE
    )
  } else if (response_scale == "percent" ||
    (response_scale == "proportion" && is.finite(max_y) && max_y > 1.5)) {
    divided_by <- 100
    scale_detected <- "percent"
    long$y <- long$y / 100
    warning(
      "Detected percent-scale responses (max > 1.5); divided y by 100 to map to [0, 1].",
      call. = FALSE
    )
  }

  # --- clamp mild out-of-range (loud, named counts) ------------------------
  n_clamped_hi <- sum(long$y > 1, na.rm = TRUE)
  n_clamped_lo <- sum(long$y < 0, na.rm = TRUE)
  if (n_clamped_hi > 0L) long$y[!is.na(long$y) & long$y > 1] <- 1
  if (n_clamped_lo > 0L) long$y[!is.na(long$y) & long$y < 0] <- 0
  if (n_clamped_hi > 0L || n_clamped_lo > 0L) {
    warning(
      sprintf(
        "Clamped out-of-range y to [0, 1]: %d value(s) > 1 set to 1, %d value(s) < 0 set to 0.",
        n_clamped_hi, n_clamped_lo
      ),
      call. = FALSE
    )
  }

  # --- error on residual NA in y -------------------------------------------
  if (anyNA(long$y)) {
    stop(
      sprintf("`%s` contains %d NA value(s) after coercion; remove or impute them.",
        y_var, sum(is.na(long$y))),
      call. = FALSE
    )
  }
  if (anyNA(long$id) || anyNA(long$x)) {
    stop("`id` and `x` must not contain NA.", call. = FALSE)
  }

  list(
    data = long,
    coercion_info = list(
      divided_by = divided_by,
      n_clamped_hi = as.integer(n_clamped_hi),
      n_clamped_lo = as.integer(n_clamped_lo),
      scale_detected = scale_detected
    )
  )
}
```

- [ ] **Step 4: Run the test, expect PASS.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-dd-validate.R")'
```
Expected output: all blocks pass; `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 7 ]` (warnings are
caught by `expect_warning`, so they do not surface as WARN).

- [ ] **Step 5: Add a NEWS note documenting the loud-coercion behavior** (spec §5 requires it
  be documented in NEWS). If `NEWS.md` does not exist, create it; otherwise prepend a 0.4.0
  bullet.

```
# NEWS.md (top of file)
# beezdiscounting 0.4.0

* New IP-family mixed-effects discounting tier (`fit_dd_tmb()`, SLT-beta and
  Gaussian families). The input validator `.dd_validate_ip()` coerces
  percent- and amount-scaled indifference points to `[0, 1]` and clamps mild
  out-of-range values, **always emitting a warning that names how many values
  were divided or clamped**. Proportion data in `[0, 1]` (including exact 0
  and 1) pass through silently.
```

- [ ] **Step 6: Commit (targeted add).**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add R/dd-validate.R tests/testthat/test-dd-validate.R NEWS.md && \
git commit -m "feat(validate): add .dd_validate_ip with loud coerce/clamp + NEWS

Coerce percent (>1.5 or response_scale='percent', /100) and amount
(/ll) responses to [0,1]; clamp mild overshoot (>1->1, <0->0) with a
warning that names the counts; error on missing id/x/y columns or NA y
after coercion. Document the loud-coercion behavior in roxygen and NEWS
per the design spec.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task F.4: `.dd_normalize_re` — random-effects formula normalizer (`k ~ 1`)

Normalize the `random_effects` argument to a canonical block list. The MVP accepts only the
single intercept-only block `k ~ 1` (natural name `"k"`; the model fits on `log k`
internally). Richer specs (slopes, multiple blocks, non-`k` LHS, pdMat objects) `stop()` with
a clear message pointing to the fast-follow. This is a simplification of beezdemand's
`.normalize_re_input` (`beezdemand/R/random-effects-utils.R:58-115`) reduced to the
single-block discounting case.

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-re-utils.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-re-parser.R`

- [ ] **Step 1: Write failing test `test-dd-tmb-re-parser.R`** (accepts `k ~ 1`; rejects
  slopes, non-k LHS, two-LHS, character vectors, pdMat-style lists).

```r
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
```

- [ ] **Step 2: Run the test, expect FAIL** (normalizer absent).

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-dd-tmb-re-parser.R")'
```
Expected output: `could not find function ".dd_normalize_re"`; failures reported.

- [ ] **Step 3: Implement `R/dd-re-utils.R`** (complete; MVP single intercept-only `k` block,
  everything richer `stop()`s).

```r
# R/dd-re-utils.R

#' Normalize the `random_effects` argument to a canonical block list
#'
#' MVP scope: only a single intercept-only random effect on the discount rate,
#' supplied as the formula `k ~ 1`, is supported. The natural parameter name is
#' `"k"`; the model is fit on `log k` internally. Anything richer (random
#' slopes, intercept suppression, multiple LHS parameters, non-`k` LHS,
#' pdMat / list / pdBlocked objects) is rejected with a message pointing to the
#' fast-follow. This is the single-random-intercept special case of
#' `beezdemand`'s `.normalize_re_input()`
#' (R/random-effects-utils.R:58-115).
#'
#' @param random_effects A two-sided formula; the only supported value is
#'   `k ~ 1`.
#' @param data Optional data frame (unused in the MVP; accepted for signature
#'   compatibility with the richer fast-follow parser).
#' @return list(source = "formula", blocks = list(list(param = "k",
#'   terms = "(Intercept)", pdmat_class = "pdDiag", formula = k ~ 1, dim = 1L))).
#' @keywords internal
#' @noRd
.dd_normalize_re <- function(random_effects, data = NULL) {
  if (!inherits(random_effects, "formula")) {
    stop(
      "`random_effects` must be a formula; the only supported value in this ",
      "release is `k ~ 1` (a single random intercept on log k).",
      call. = FALSE
    )
  }
  if (length(random_effects) != 3L) {
    stop(
      "`random_effects` must be a two-sided formula `k ~ 1`.",
      call. = FALSE
    )
  }

  lhs_vars <- all.vars(random_effects[[2]])
  if (!identical(lhs_vars, "k")) {
    stop(
      "The random-effects LHS must be exactly `k` (got ",
      paste(shQuote(lhs_vars), collapse = ", "),
      "). Subject-random phi and 2-parameter random effects are out of scope ",
      "in this release.",
      call. = FALSE
    )
  }

  rhs <- random_effects[[3]]
  tt <- stats::terms(stats::as.formula(paste("~", deparse1(rhs))))
  has_intercept <- attr(tt, "intercept") == 1L
  n_terms <- length(attr(tt, "term.labels"))
  if (!has_intercept || n_terms > 0L) {
    stop(
      "Only a single intercept-only random effect `k ~ 1` is supported in ",
      "this release; random slopes are a fast-follow.",
      call. = FALSE
    )
  }

  list(
    source = "formula",
    blocks = list(list(
      param = "k",
      terms = "(Intercept)",
      pdmat_class = "pdDiag",
      formula = k ~ 1,
      dim = 1L
    ))
  )
}
```

- [ ] **Step 4: Run the test, expect PASS.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'testthat::test_file("tests/testthat/test-dd-tmb-re-parser.R")'
```
Expected output: all blocks pass; `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 7 ]`.

- [ ] **Step 5: Run the full pure-R test slice (lint + all four new test files) to confirm the
  foundation is green together.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
Rscript -e 'for (f in c("tests/testthat/test-sltb-density.R","tests/testthat/test-dd-param-space.R","tests/testthat/test-dd-validate.R","tests/testthat/test-dd-tmb-re-parser.R")) testthat::test_file(f)'
```
Expected output: every file reports `FAIL 0`; combined PASS count = 35.

- [ ] **Step 6: Commit (targeted add).**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add R/dd-re-utils.R tests/testthat/test-dd-tmb-re-parser.R && \
git commit -m "feat(re): add .dd_normalize_re for the k ~ 1 single random intercept

Accept only the intercept-only formula k ~ 1 (natural name k, fit on
log k internally) and return the canonical single-block representation;
reject slopes, intercept suppression, non-k LHS, two-parameter LHS, and
non-formula input with messages pointing to the fast-follows. Simplified
from beezdemand .normalize_re_input.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
## Phase T — TMB core (C++ template + build)

This phase builds the compiled heart of the package: the `MixedDiscounting` TMB template
(`eqn_type × family` switch with the SLT-beta and Gaussian densities), its dispatcher with
`TMB_LIB_INIT R_init_beezdiscounting`, the `src/Makevars*` build files, the `DESCRIPTION`
edits that turn `beezdiscounting` into a compiled package, and the compile gate that
cross-checks the C++ `-nll` against the R SLT log-density to `1e-8`.

**Prerequisite assumed from an earlier phase:** `R/dd-param-space.R` defines the pure-R SLT
log-density helper `.dd_slt_logpdf(y, mu, phi, s = 1.0000001, l = 1e-8)` (the verbatim port of
`slt_logpdf` from `dev/sltb-verification/verify_sltb.R`). The cross-check test in Task T.5
**uses** that helper as the source of truth. If it is not yet present when this phase runs,
Task T.5 includes the minimal inline definition needed for the gate (see Step 1 of T.5), but
the canonical home is `R/dd-param-space.R`.

---

### Task T.1: DESCRIPTION + `.Rbuildignore` — make `beezdiscounting` a compiled package

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/DESCRIPTION` (Imports block ends at the line `tidyr`; Suggests block is the single line `testthat (>= 3.0.0)`)
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/.Rbuildignore` (already has `^src/.*\.o$` and `^src/.*\.so$`; missing `.dll` and `symbols.rds`)

- [ ] **Step 1: Add `LinkingTo`, extend `Imports`, extend `Suggests`, add build fields to `DESCRIPTION`.**

  Edit the `Imports:` block to append three packages (keep existing entries; alphabetical not required but append cleanly):

  ```
  Imports:
      beezdemand,
      broom,
      dplyr,
      emmeans,
      ggplot2,
      gtools,
      magrittr,
      minpack.lm,
      psych,
      purrr,
      RcppEigen,
      stringr,
      tibble,
      tidyr,
      TMB (>= 1.9.0)
  ```

  Replace the `Suggests:` block with:

  ```
  Suggests:
      knitr,
      rmarkdown,
      testthat (>= 3.0.0)
  ```

  Add these fields immediately after the `Suggests:` block (before `LazyData: true` or at end of file):

  ```
  LinkingTo:
      RcppEigen,
      TMB
  VignetteBuilder: knitr
  SystemRequirements: GNU make
  Config/testthat/edition: 3
  ```

  Bump the version line to:

  ```
  Version: 0.4.0
  ```

- [ ] **Step 2: Verify/extend `.Rbuildignore` for compiled artifacts.**

  The file already contains `^src/.*\.o$` and `^src/.*\.so$`. Add the two missing artifact
  patterns so Windows DLLs and the TMB symbols table are not shipped. After the existing
  `^src/.*\.so$` line, add:

  ```
  ^src/.*\.dll$
  ^src/symbols\.rds$
  ```

- [ ] **Step 3: Run — confirm DESCRIPTION parses and fields are present.**

  Run:
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'd <- read.dcf("DESCRIPTION"); stopifnot("LinkingTo" %in% colnames(d)); stopifnot(grepl("TMB", d[,"LinkingTo"])); stopifnot(grepl("TMB", d[,"Imports"])); stopifnot(grepl("emmeans", d[,"Imports"])); stopifnot(grepl("RcppEigen", d[,"Imports"])); stopifnot(d[,"Version"] == "0.4.0"); stopifnot(grepl("GNU make", d[,"SystemRequirements"])); cat("DESCRIPTION OK\n")'
  ```
  Expected output:
  ```
  DESCRIPTION OK
  ```

- [ ] **Step 4: Commit.**
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting
  git add DESCRIPTION .Rbuildignore
  git commit -m "build(tmb): add TMB/RcppEigen/emmeans deps and compiled-pkg fields

Add LinkingTo: RcppEigen, TMB; Imports: TMB (>= 1.9.0), RcppEigen, emmeans;
Suggests: knitr, rmarkdown; VignetteBuilder: knitr; SystemRequirements: GNU make.
Bump to 0.4.0. Ignore src .dll and symbols.rds build artifacts.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task T.2: `src/Makevars` and `src/Makevars.win`

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/src/Makevars`
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/src/Makevars.win`

beezdemand ships an **empty** `Makevars`/`Makevars.win` (verified: both files exist but contain
no `PKG_*` lines — TMB/RcppEigen flags come entirely from `LinkingTo` headers). We mirror that:
the files exist (so R's build system treats `src/` as a compiled dir and the platform `Makevars`
is selected on Windows) but carry only a documenting comment. No `PKG_CPPFLAGS` is needed
because TMB and RcppEigen inject their include paths through `LinkingTo`.

- [ ] **Step 1: Create `src/Makevars`.**

  ```make
  # Compiler/linker flags for the beezdiscounting TMB template.
  # TMB and RcppEigen supply their include paths via LinkingTo in DESCRIPTION,
  # so no PKG_CPPFLAGS/PKG_LIBS are required here (mirrors beezdemand 0.3.0).
  ```

- [ ] **Step 2: Create `src/Makevars.win`.**

  ```make
  # Windows compiler/linker flags for the beezdiscounting TMB template.
  # TMB and RcppEigen supply their include paths via LinkingTo in DESCRIPTION,
  # so no PKG_CPPFLAGS/PKG_LIBS are required here (mirrors beezdemand 0.3.0).
  ```

- [ ] **Step 3: Run — confirm files exist and are non-binary.**

  Run:
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting && test -f src/Makevars && test -f src/Makevars.win && file src/Makevars src/Makevars.win && echo "MAKEVARS OK"
  ```
  Expected output (last line):
  ```
  MAKEVARS OK
  ```

- [ ] **Step 4: Commit.**
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting
  git add src/Makevars src/Makevars.win
  git commit -m "build(tmb): add empty src/Makevars and src/Makevars.win

LinkingTo (TMB, RcppEigen) supplies all include paths; the platform Makevars
files exist only to mark src/ as a compiled directory (mirrors beezdemand).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task T.3: `src/beezdiscounting_TMBExports.cpp` — DLL init + `MixedDiscounting` dispatch

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/src/beezdiscounting_TMBExports.cpp`

This is a near-verbatim port of `/Users/brent/Dropbox/GIT/beezdemand/src/beezdemand_TMBExports.cpp`
(lines 1–37) with three adaptations: (a) `TMB_LIB_INIT R_init_beezdemand` →
`TMB_LIB_INIT R_init_beezdiscounting`; (b) the include list is reduced to the single
`MixedDiscounting.h`; (c) the dispatch branches collapse to one `if (model == "MixedDiscounting")`.

- [ ] **Step 1: Create the dispatcher.**

  ```cpp
  // beezdiscounting TMB model dispatcher
  // Registers the MixedDiscounting TMB model for the beezdiscounting package.

  #define TMB_LIB_INIT R_init_beezdiscounting
  #include <TMB.hpp>
  #include "MixedDiscounting.h"

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_STRING(model);
    if (model == "MixedDiscounting") {
      return MixedDiscounting(this);
    } else {
      error("Unknown model");
    }
    return Type(0);
  }
  ```

  (No standalone run/commit here — `MixedDiscounting.h` does not yet exist, so the package will
  not compile until Task T.4. Commit both together at the end of T.4's compile gate.)

---

### Task T.4: `src/MixedDiscounting.h` — the full C++ template (eqn_type × family)

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/src/MixedDiscounting.h`

This is **novel** C++ (the SLT-beta density has no analog in beezdemand). Structure follows the
single-RE pattern of `HurdleDemand2RE.h` (standardized `u`, non-centered prior) and the
`eqn_type` runtime switch / `CppAD::CondExp` guards of `MixedDemand.h`. Data/parameter names
match the PINNED CONTRACT R-side lists exactly. The SLT density is the verified math from
`dev/sltb-verification/verify_sltb.R::slt_logpdf` (`s = 1.0000001`, `l = 1e-8`), ported to
`Type` with `pbeta` atomic.

Key contract facts encoded here:
- `parameters = list(beta_k=<ncol(X)>, log_sigma_u=<scalar>, log_aux=<scalar>, u=<n_subjects x 1>)`.
- `log_aux` is the single generic auxiliary scalar: `phi = exp(log_aux)` when `family==0` (sltb);
  `sigma_e = exp(log_aux)` when `family==1` (gaussian).
- linear predictor: `log_k_i = (X.row(i)*beta_k).sum() + sigma_u * u(subj,0)`; `k_i = exp(log_k_i)`.
- `mu_i`: mazur (`eqn_type==0`) `1/(1+k_i*x_i)`; exponential (`eqn_type==1`) `exp(-k_i*x_i)`;
  guarded to `[1e-6, 1-1e-6]` via `CppAD::CondExpLt/Gt`.
- RE prior: `nll -= dnorm(u(i,0), 0, 1, true)` (non-centered; `sigma_u` multiplies in predictor).
- ADREPORT: `beta_k`, `sigma_u = exp(log_sigma_u)`, and `aux` (`phi` or `sigma_e`).

- [ ] **Step 1: Create the full template.**

  ```cpp
  /// @file MixedDiscounting.h
  /// IP-family Mixed-Effects Discounting via TMB (beezdiscounting 0.4.0)
  /// =============================================================================
  ///
  /// MODEL STRUCTURE:
  /// One-parameter discounting (Mazur hyperbolic or exponential) for indifference
  /// proportions y in [0,1], with a random intercept on log k and between-subject
  /// fixed effects via design matrix X. Two observation families dispatch at
  /// runtime: scale-location-truncated beta (SLT-beta) and Gaussian.
  ///
  ///   eqn_type: 0 = mazur        mu = 1 / (1 + k*x)
  ///             1 = exponential  mu = exp(-k*x)
  ///   family:   0 = sltb         y ~ SLTBeta(mu, phi)      aux = phi
  ///             1 = gaussian     y ~ N(mu, sigma_e)        aux = sigma_e
  ///
  /// LINEAR PREDICTOR (log k scale, single random intercept, non-centered):
  ///   log_k_i = (X.row(i) * beta_k).sum() + sigma_u * u(subj, 0)
  ///   k_i     = exp(log_k_i)
  ///   sigma_u = exp(log_sigma_u),  u(i,0) ~ N(0,1)
  ///
  /// AUXILIARY: log_aux is a single generic scalar.
  ///   family == 0 (sltb):     phi     = exp(log_aux)
  ///   family == 1 (gaussian): sigma_e = exp(log_aux)
  ///
  /// SLT-BETA DENSITY (constants s = 1.0000001, l = 1e-8; verified in
  /// dev/sltb-verification/verify_sltb.R):
  ///   a = mu*phi, b = (1-mu)*phi
  ///   log f = lgamma(a+b) - lgamma(a) - lgamma(b)
  ///         + (a-1)*log(y/s + l) + (b-1)*log(1 - (y/s + l))
  ///         - log(s)
  ///         - log( pbeta(1/s + l, a, b) - pbeta(l, a, b) )
  /// The truncation normalizer Z = pbeta(1/s+l,a,b) - pbeta(l,a,b) is load-bearing.
  ///
  /// MU is guarded to [1e-6, 1-1e-6] via CppAD::CondExp (exponential underflows to
  /// 0 at long delays). subject_id is 0-indexed from R.
  /// =============================================================================

  #undef TMB_OBJECTIVE_PTR
  #define TMB_OBJECTIVE_PTR obj

  template <class Type>
  Type MixedDiscounting(objective_function<Type>* obj) {
    // =========================================================================
    // DATA
    // =========================================================================
    DATA_VECTOR(y);               // indifference proportions in [0,1], length n_obs
    DATA_VECTOR(x);               // delays, length n_obs
    DATA_IVECTOR(subject_id);     // 0-indexed subject index per observation
    DATA_MATRIX(X);               // n_obs x ncol(X) fixed-effect design for log k
    DATA_INTEGER(eqn_type);       // 0 = mazur, 1 = exponential
    DATA_INTEGER(family);         // 0 = sltb, 1 = gaussian
    DATA_INTEGER(n_obs);
    DATA_INTEGER(n_subjects);

    // =========================================================================
    // PARAMETERS
    // =========================================================================
    PARAMETER_VECTOR(beta_k);     // length ncol(X): fixed effects on log k
    PARAMETER(log_sigma_u);       // log SD of random intercept on log k
    PARAMETER(log_aux);           // log phi (sltb) OR log sigma_e (gaussian)
    PARAMETER_MATRIX(u);          // n_subjects x 1 standardized random intercepts

    // =========================================================================
    // TRANSFORM SCALARS
    // =========================================================================
    Type sigma_u = exp(log_sigma_u);
    Type aux = exp(log_aux);      // phi (family 0) or sigma_e (family 1)

    // SLT-beta fixed constants (verified reference values).
    Type s = Type(1.0000001);
    Type l = Type(1e-8);

    // mu guard bounds.
    Type mu_lo = Type(1e-6);
    Type mu_hi = Type(1.0) - Type(1e-6);

    Type log2pi = log(Type(2.0) * M_PI);

    // =========================================================================
    // NEGATIVE LOG-LIKELIHOOD
    // =========================================================================
    Type nll = Type(0.0);

    // Standard-normal prior on the standardized random intercepts (non-centered;
    // sigma_u enters multiplicatively in the linear predictor below).
    for (int i = 0; i < n_subjects; i++) {
      nll -= dnorm(u(i, 0), Type(0.0), Type(1.0), true);
    }

    // Data likelihood.
    for (int i = 0; i < n_obs; i++) {
      int subj = subject_id(i);

      // Linear predictor on log k: fixed effects + random intercept.
      vector<Type> x_i = X.row(i);
      Type log_k_i = (x_i * beta_k).sum() + sigma_u * u(subj, 0);
      Type k_i = exp(log_k_i);

      // Mean via the discounting function (identity link).
      Type mu_raw;
      if (eqn_type == 0) {
        // Mazur hyperbolic.
        mu_raw = Type(1.0) / (Type(1.0) + k_i * x(i));
      } else {
        // Exponential.
        mu_raw = exp(-k_i * x(i));
      }

      // Guard mu to [mu_lo, mu_hi] (exponential underflow at long delays).
      // Never branch on a Type-valued condition with `if`; use CppAD::CondExp.
      Type mu = CppAD::CondExpLt(mu_raw, mu_lo, mu_lo, mu_raw);
      mu = CppAD::CondExpGt(mu, mu_hi, mu_hi, mu);

      if (family == 0) {
        // ------------------------------------------------------------------
        // SLT-beta density. aux = phi.
        // ------------------------------------------------------------------
        Type phi = aux;
        Type a = mu * phi;
        Type b = (Type(1.0) - mu) * phi;

        Type yt = y(i) / s + l;                  // scaled-location transform
        Type Z = pbeta(Type(1.0) / s + l, a, b)  // truncation normalizer
                 - pbeta(l, a, b);

        Type logf = lgamma(a + b) - lgamma(a) - lgamma(b)
                    + (a - Type(1.0)) * log(yt)
                    + (b - Type(1.0)) * log(Type(1.0) - yt)
                    - log(s)
                    - log(Z);
        nll -= logf;

      } else {
        // ------------------------------------------------------------------
        // Gaussian density. aux = sigma_e.
        // ------------------------------------------------------------------
        Type sigma_e = aux;
        Type resid = (y(i) - mu) / sigma_e;
        nll -= -log_aux - Type(0.5) * log2pi - Type(0.5) * resid * resid;
      }
    }

    // =========================================================================
    // ADREPORT
    // =========================================================================
    ADREPORT(beta_k);     // fixed effects on log k
    ADREPORT(sigma_u);    // SD of random intercept on log k
    ADREPORT(aux);        // phi (family 0) or sigma_e (family 1)

    return nll;
  }

  #undef TMB_OBJECTIVE_PTR
  #define TMB_OBJECTIVE_PTR this
  ```

  Note on the Gaussian normalizing constant: `-log_aux` equals `-log(sigma_e)` (since
  `sigma_e = exp(log_aux)`), matching beezdemand's `-logsigma_e - 0.5*log2pi - 0.5*resid^2`
  pattern (`MixedDemand.h:241`).

---

### Task T.5: Compile gate + C++ `-nll` == R SLT log-density cross-check (to 1e-8)

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-tmb-compile-gate.R`
- (No new R; relies on `R/dd-param-space.R::.dd_slt_logpdf` if present, else inlines it in the test — see Step 1.)

This is the standard TMB compile gate: compile + link the template once, build a `MakeADFun`
object for a no-RE / fixed-effect-only configuration, and assert that `obj$fn(par)` (the C++
`-nll`) equals `sum(-.dd_slt_logpdf(y, mu, phi))` (the verified R density) to `1e-8`. We test
**both families**: SLTB (against `.dd_slt_logpdf`) and Gaussian (against `dnorm(..., log=TRUE)`).

To isolate the **observation density** (the load-bearing novel math) from the RE machinery, the
gate fixes `u = 0` and `sigma_u` arbitrary (RE drops out: `sigma_u * 0 = 0`) and uses a single
subject with an intercept-only `X` (so `log_k = beta_k[1]`, i.e. `k = exp(beta_k[1])`). With
`u(0,0)=0` the RE prior term is `-dnorm(0,0,1,log=TRUE) = 0.5*log(2*pi)` per subject — the test
accounts for this constant explicitly so the comparison is exact.

- [ ] **Step 1: Write the failing test (template not yet compiled).**

  ```r
  # tests/testthat/test-tmb-compile-gate.R
  # Compile gate: the C++ MixedDiscounting -nll must equal the verified R SLT
  # log-density (and the Gaussian density) to 1e-8, isolating the observation
  # term by zeroing the random effects.

  # Verified R SLT log-density (source of truth: dev/sltb-verification/verify_sltb.R
  # slt_logpdf). Canonical home is R/dd-param-space.R::.dd_slt_logpdf; defined
  # locally here so the gate is self-contained if param-space lands later.
  .gate_slt_logpdf <- function(y, mu, phi, s = 1.0000001, l = 1e-8) {
    a <- mu * phi
    b <- (1 - mu) * phi
    lgamma(a + b) - lgamma(a) - lgamma(b) +
      (a - 1) * log(y / s + l) + (b - 1) * log(1 - (y / s + l)) -
      log(s) - log(pbeta(1 / s + l, a, b) - pbeta(l, a, b))
  }

  # Build a single-subject, intercept-only TMB data/parameter set and the
  # matching R-side prediction for a given equation/family.
  .gate_setup <- function(eqn_type, family, beta0, aux) {
    delays <- c(7, 30, 180, 365, 730, 1460, 2920)
    y <- c(0.95, 0.80, 0.50, 0.35, 0.20, 0.08, 0.00)
    n_obs <- length(delays)
    k <- exp(beta0)
    mu_raw <- if (eqn_type == 0L) 1 / (1 + k * delays) else exp(-k * delays)
    mu <- pmin(pmax(mu_raw, 1e-6), 1 - 1e-6)
    list(
      data = list(
        model = "MixedDiscounting",
        y = y, x = delays,
        subject_id = rep(0L, n_obs),
        X = matrix(1, nrow = n_obs, ncol = 1L),
        eqn_type = eqn_type, family = family,
        n_obs = n_obs, n_subjects = 1L
      ),
      parameters = list(
        beta_k = beta0,
        log_sigma_u = log(0.5),     # arbitrary; RE drops out because u = 0
        log_aux = log(aux),
        u = matrix(0, nrow = 1L, ncol = 1L)
      ),
      mu = mu, y = y
    )
  }

  describe("MixedDiscounting compile gate", {
    skip_on_cran()
    skip_if_not_installed("TMB")

    it("compiles and links the MixedDiscounting template", {
      # Loading the installed package DLL is enough; if devtools::load_all was
      # used, the DLL is already loaded. getLoadedDLLs() must list beezdiscounting.
      expect_true("beezdiscounting" %in% names(getLoadedDLLs()))
    })

    it("C++ -nll matches the R SLT log-density to 1e-8 (mazur, sltb)", {
      g <- .gate_setup(eqn_type = 0L, family = 0L, beta0 = log(0.01), aux = 8)
      obj <- TMB::MakeADFun(data = g$data, parameters = g$parameters,
                            random = "u", DLL = "beezdiscounting", silent = TRUE)
      par <- obj$par
      cpp_nll <- obj$fn(par)
      re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))           # one subject, u = 0
      r_nll <- sum(-.gate_slt_logpdf(g$y, g$mu, 8)) + re_prior
      expect_equal(as.numeric(cpp_nll), r_nll, tolerance = 1e-8)
    })

    it("C++ -nll matches the R SLT log-density to 1e-8 (exponential, sltb)", {
      g <- .gate_setup(eqn_type = 1L, family = 0L, beta0 = log(2e-4), aux = 8)
      obj <- TMB::MakeADFun(data = g$data, parameters = g$parameters,
                            random = "u", DLL = "beezdiscounting", silent = TRUE)
      cpp_nll <- obj$fn(obj$par)
      re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
      r_nll <- sum(-.gate_slt_logpdf(g$y, g$mu, 8)) + re_prior
      expect_equal(as.numeric(cpp_nll), r_nll, tolerance = 1e-8)
    })

    it("C++ -nll matches the R Gaussian density to 1e-8 (mazur, gaussian)", {
      g <- .gate_setup(eqn_type = 0L, family = 1L, beta0 = log(0.01), aux = 0.1)
      obj <- TMB::MakeADFun(data = g$data, parameters = g$parameters,
                            random = "u", DLL = "beezdiscounting", silent = TRUE)
      cpp_nll <- obj$fn(obj$par)
      re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
      r_nll <- sum(-dnorm(g$y, g$mu, 0.1, log = TRUE)) + re_prior
      expect_equal(as.numeric(cpp_nll), r_nll, tolerance = 1e-8)
    })
  })
  ```

  Run (expect FAIL — template not compiled, DLL not loaded):
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-tmb-compile-gate.R")'
  ```
  Expected output (the `src/*.cpp` will fail to compile until T.3+T.4 land, or the DLL won't
  be loaded): a compile error from `devtools::load_all` OR a failing first `it()` with
  `"beezdiscounting" %in% names(getLoadedDLLs())` not TRUE. Either way: NOT all green.

- [ ] **Step 2: Implement — ensure T.3 (`beezdiscounting_TMBExports.cpp`) and T.4 (`MixedDiscounting.h`) are in place, then compile.**

  Both source files were authored in Tasks T.3 and T.4. Compile the package's `src/`:
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'TMB::compile("src/beezdiscounting_TMBExports.cpp"); cat("COMPILE OK\n")'
  ```
  Expected output (after compiler chatter):
  ```
  COMPILE OK
  ```
  If `pbeta` / `lgamma` raise an unknown-symbol error, confirm `#include <TMB.hpp>` precedes
  `#include "MixedDiscounting.h"` in `beezdiscounting_TMBExports.cpp` (it does) — TMB's atomic
  `pbeta` and `lgamma` come from `TMB.hpp`.

- [ ] **Step 3: Run the gate via `devtools::load_all` (which compiles + loads the DLL), expect PASS.**

  Run:
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-tmb-compile-gate.R")'
  ```
  Expected output:
  ```
  [ FAIL 0 | WARN 0 | SKIP 0 | PASS 4 ]
  ```
  (4 PASS = compile/DLL-loaded + mazur-sltb + exponential-sltb + mazur-gaussian. SKIP 0 because
  TMB is installed locally; on CRAN the `skip_on_cran()` / `skip_if_not_installed("TMB")` guards
  fire and these become SKIPs.)

- [ ] **Step 4: Commit the template, dispatcher, and gate together.**
  ```
  cd /Users/brent/Dropbox/GIT/beezdiscounting
  git add src/MixedDiscounting.h src/beezdiscounting_TMBExports.cpp tests/testthat/test-tmb-compile-gate.R
  git commit -m "feat(tmb): MixedDiscounting C++ template + compile gate

Add src/MixedDiscounting.h (eqn_type x family: mazur/exponential, SLT-beta/
gaussian), the beezdiscounting_TMBExports.cpp dispatcher with
TMB_LIB_INIT R_init_beezdiscounting, and a compile gate asserting the C++ -nll
equals the verified R SLT log-density (and Gaussian density) to 1e-8.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---
## Phase P — Fit pipeline & fit_dd_tmb()

This phase wires the `.dd_tmb_*` helpers and the public `fit_dd_tmb()`. It depends on
earlier phases having created: `R/dd-validate.R` (`.dd_validate_ip`), `R/dd-tmb.R`
(`.dd_normalize_re`), and the compiled DLL `"beezdiscounting"` exposing the
`"MixedDiscounting"` model (from the C++/TMB phase). All compiled-template tests start
with `skip_on_cran()` and `skip_if_not_installed("TMB")`. Tests are testthat 3e, BDD
`describe()/it()`. Use TARGETED `git add` only.

### Task P.1: `.dd_tmb_prepare_data()` — 0-indexed subject map, NA drop, row-order coherence

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append helper)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (new file; this task adds the prepare-data `describe` block)

- [ ] **Step 1: Write failing test for `.dd_tmb_prepare_data`.**

Create `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R`:

```r
describe(".dd_tmb_prepare_data()", {
  it("0-indexes subject_id and aligns vectors with subject_levels", {
    dat <- data.frame(
      id = factor(c("b", "b", "a", "a"), levels = c("b", "a")),
      x = c(7, 30, 7, 30),
      y = c(0.9, 0.5, 0.8, 0.4)
    )
    prep <- .dd_tmb_prepare_data(dat, y_var = "y", x_var = "x", id_var = "id")
    expect_equal(prep$subject_levels, c("a", "b"))
    expect_equal(prep$n_subjects, 2L)
    expect_equal(prep$n_obs, 4L)
    expect_equal(min(prep$subject_id), 0L)
    expect_equal(max(prep$subject_id), 1L)
    # subject_id maps each row to the sorted-level index
    expect_equal(prep$subject_id, c(1L, 1L, 0L, 0L))
    expect_type(prep$subject_id, "integer")
  })

  it("drops rows with NA in y/x/id and keeps x,y,subject_id row-coherent", {
    dat <- data.frame(
      id = c("a", "a", "b", "b"),
      x = c(7, 30, 7, NA),
      y = c(0.9, NA, 0.8, 0.4)
    )
    prep <- .dd_tmb_prepare_data(dat, y_var = "y", x_var = "x", id_var = "id")
    expect_equal(prep$n_obs, 2L)
    expect_equal(prep$y, c(0.9, 0.8))
    expect_equal(prep$x, c(7, 7))
    expect_equal(nrow(prep$data), 2L)
    # surviving rows keep their parallel arrays aligned
    expect_equal(length(prep$x), length(prep$subject_id))
  })

  it("renames caller columns to canonical id/x/y in $data", {
    dat <- data.frame(pid = "p1", delay = 7, ip = 0.5)
    prep <- .dd_tmb_prepare_data(dat, y_var = "ip", x_var = "delay", id_var = "pid")
    expect_named(prep$data, c("id", "x", "y"))
  })
})
```

- [ ] **Step 2: Run the test — expect failure (`.dd_tmb_prepare_data` not found).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R")'
```
Expected output: error `could not find function ".dd_tmb_prepare_data"` (all `describe` blocks fail).

- [ ] **Step 3: Implement `.dd_tmb_prepare_data` in `R/dd-tmb.R`.**

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
#' Prepare data for the TMB mixed-effects discounting model
#'
#' Maps caller column names to canonical `id`/`x`/`y`, drops rows with NA in any
#' modeling column, and builds a 0-indexed `subject_id` aligned to
#' `subject_levels` (the C++ template indexes `u(subject_id, 0)` from 0).
#'
#' @param data Data frame already coerced/clamped by [.dd_validate_ip()].
#' @param y_var,x_var,id_var Character column names.
#' @return A list with `y`, `x`, `subject_id` (0-indexed integer),
#'   `subject_levels`, `n_subjects`, `n_obs`, and `data` (cleaned long df with
#'   canonical columns `id`, `x`, `y`).
#' @keywords internal
.dd_tmb_prepare_data <- function(data, y_var = "y", x_var = "x", id_var = "id") {
  ids <- data[[id_var]]
  x <- as.numeric(data[[x_var]])
  y <- as.numeric(data[[y_var]])

  keep <- !is.na(ids) & !is.na(x) & !is.na(y)
  ids <- ids[keep]
  x <- x[keep]
  y <- y[keep]

  if (length(y) == 0L) {
    stop("No complete cases remain after dropping NA rows.", call. = FALSE)
  }

  # 0-indexed subject map, levels sorted for reproducibility
  subject_levels <- sort(unique(as.character(ids)))
  n_subjects <- length(subject_levels)
  subject_map <- stats::setNames(seq_along(subject_levels) - 1L, subject_levels)
  subject_id <- as.integer(subject_map[as.character(ids)])

  cleaned <- data.frame(
    id = as.character(ids),
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )

  list(
    y = y,
    x = x,
    subject_id = subject_id,
    subject_levels = subject_levels,
    n_subjects = n_subjects,
    n_obs = length(y),
    data = cleaned
  )
}
```

- [ ] **Step 4: Run the test — expect pass for the prepare-data block.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R")'
```
Expected output: the `.dd_tmb_prepare_data()` `describe` block passes (3 `it` green); later blocks still fail (helpers not yet written).

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): add .dd_tmb_prepare_data (0-indexed subject map, NA drop)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task P.2: `.dd_tmb_build_design()` — model.matrix for log-k fixed effects

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append helper)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (add `describe` block)

- [ ] **Step 1: Write failing test for `.dd_tmb_build_design`.**

Append to `tests/testthat/test-fit_dd_tmb.R`:

```r
describe(".dd_tmb_build_design()", {
  it("returns an intercept-only X when no factors/covariates given", {
    dat <- data.frame(id = "a", x = 7, y = 0.5)
    d <- .dd_tmb_build_design(dat)
    expect_equal(colnames(d$X), "(Intercept)")
    expect_equal(ncol(d$X), 1L)
    expect_equal(d$rhs, "~ 1")
  })

  it("expands a single between-subject factor into contrast columns", {
    dat <- data.frame(
      id = rep(c("a", "b", "c", "d"), each = 1),
      x = 7,
      y = 0.5,
      grp = factor(c("ctrl", "ctrl", "trt", "trt"))
    )
    d <- .dd_tmb_build_design(dat, factors = "grp")
    expect_equal(ncol(d$X), 2L)
    expect_true("grptrt" %in% colnames(d$X))
    expect_match(d$rhs, "grp")
  })

  it("adds an interaction term when factor_interaction = TRUE", {
    dat <- data.frame(
      id = letters[1:4], x = 7, y = 0.5,
      a = factor(c("x", "x", "z", "z")),
      b = factor(c("p", "q", "p", "q"))
    )
    d <- .dd_tmb_build_design(dat, factors = c("a", "b"),
                              factor_interaction = TRUE)
    expect_true(any(grepl(":", colnames(d$X))))
  })

  it("adds continuous covariate main-effect columns", {
    dat <- data.frame(id = letters[1:3], x = 7, y = 0.5, age = c(20, 30, 40))
    d <- .dd_tmb_build_design(dat, continuous_covariates = "age")
    expect_true("age" %in% colnames(d$X))
  })
})
```

- [ ] **Step 2: Run the test — expect failure (`.dd_tmb_build_design` not found).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_build_design()` block errors (function not found).

- [ ] **Step 3: Implement `.dd_tmb_build_design` in `R/dd-tmb.R`.**

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
#' Build the fixed-effect design matrix for log k
#'
#' Constructs the `model.matrix` for the `log k = Xbeta` linear predictor from
#' between-subject factors, an optional pairwise interaction, and continuous
#' covariates. With no inputs the design is intercept-only.
#'
#' @param data Cleaned long data frame (canonical `id`/`x`/`y` plus any factor
#'   or covariate columns).
#' @param factors Character vector of factor column names, or `NULL`.
#' @param factor_interaction Logical; if `TRUE` and >= 2 factors, include their
#'   interaction (uses `*`); otherwise main effects (`+`).
#' @param continuous_covariates Character vector of covariate names, or `NULL`.
#' @return A list with `X` (the model matrix for log-k FE) and `rhs` (the RHS
#'   formula string used to build it).
#' @keywords internal
.dd_tmb_build_design <- function(data, factors = NULL,
                                 factor_interaction = FALSE,
                                 continuous_covariates = NULL) {
  terms <- character(0)

  if (!is.null(factors) && length(factors) > 0L) {
    if (isTRUE(factor_interaction) && length(factors) >= 2L) {
      terms <- c(terms, paste(factors, collapse = " * "))
    } else {
      terms <- c(terms, factors)
    }
  }
  if (!is.null(continuous_covariates) && length(continuous_covariates) > 0L) {
    terms <- c(terms, continuous_covariates)
  }

  rhs <- if (length(terms) == 0L) "~ 1" else paste("~", paste(terms, collapse = " + "))
  X <- stats::model.matrix(stats::as.formula(rhs), data = data)

  list(X = X, rhs = rhs)
}
```

- [ ] **Step 4: Run the test — expect pass.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_build_design()` block passes (4 `it` green).

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): add .dd_tmb_build_design (model.matrix for log-k FE)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task P.3: `.dd_tmb_build_tmb_data()` + `.dd_tmb_default_starts()`

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append two helpers)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (add `describe` block)

- [ ] **Step 1: Write failing test for the TMB data list + default starts.**

Append to `tests/testthat/test-fit_dd_tmb.R`:

```r
describe(".dd_tmb_build_tmb_data()", {
  it("assembles the TMB data list with the contract names and enum ints", {
    prep <- list(
      y = c(0.9, 0.5), x = c(7, 30), subject_id = c(0L, 0L),
      subject_levels = "a", n_subjects = 1L, n_obs = 2L,
      data = data.frame(id = "a", x = c(7, 30), y = c(0.9, 0.5))
    )
    design <- list(X = matrix(1, nrow = 2, ncol = 1,
                              dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    td <- .dd_tmb_build_tmb_data(prep, design, equation = "mazur",
                                 family = "sltb")
    expect_equal(td$model, "MixedDiscounting")
    expect_equal(td$eqn_type, 0L)   # mazur
    expect_equal(td$family, 0L)     # sltb
    expect_equal(td$n_obs, 2L)
    expect_equal(td$n_subjects, 1L)
    expect_equal(td$subject_id, c(0L, 0L))
    expect_true(is.matrix(td$X))
  })

  it("maps exponential -> 1 and gaussian -> 1", {
    prep <- list(y = 0.5, x = 7, subject_id = 0L, subject_levels = "a",
                 n_subjects = 1L, n_obs = 1L,
                 data = data.frame(id = "a", x = 7, y = 0.5))
    design <- list(X = matrix(1, 1, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    td <- .dd_tmb_build_tmb_data(prep, design, equation = "exponential",
                                 family = "gaussian")
    expect_equal(td$eqn_type, 1L)
    expect_equal(td$family, 1L)
  })
})

describe(".dd_tmb_default_starts()", {
  it("derives log_k intercept from median y at min delay (mazur)", {
    # median y at min delay = 0.8 -> mu=0.8 -> k = (1/mu - 1)/x_min
    prep <- list(
      y = c(0.8, 0.8, 0.3, 0.3), x = c(7, 7, 365, 365),
      subject_id = c(0L, 1L, 0L, 1L), subject_levels = c("a", "b"),
      n_subjects = 2L, n_obs = 4L,
      data = data.frame(id = c("a", "b", "a", "b"),
                        x = c(7, 7, 365, 365), y = c(0.8, 0.8, 0.3, 0.3))
    )
    design <- list(X = matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb")
    expect_length(st$beta_k, 1L)
    k_implied <- (1 / 0.8 - 1) / 7
    expect_equal(st$beta_k[1], log(k_implied), tolerance = 1e-8)
    expect_equal(st$log_sigma_u, log(0.5))
    expect_equal(st$log_aux, log(8))          # sltb
    expect_equal(dim(st$u), c(2L, 1L))
  })

  it("uses log(0.1) for log_aux under gaussian", {
    prep <- list(y = c(0.8, 0.3), x = c(7, 365), subject_id = c(0L, 0L),
                 subject_levels = "a", n_subjects = 1L, n_obs = 2L,
                 data = data.frame(id = "a", x = c(7, 365), y = c(0.8, 0.3)))
    design <- list(X = matrix(1, 2, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "gaussian")
    expect_equal(st$log_aux, log(0.1))
  })

  it("zero-pads beta_k for multi-column designs", {
    prep <- list(y = c(0.8, 0.5), x = c(7, 7), subject_id = c(0L, 1L),
                 subject_levels = c("a", "b"), n_subjects = 2L, n_obs = 2L,
                 data = data.frame(id = c("a", "b"), x = 7, y = c(0.8, 0.5)))
    design <- list(
      X = matrix(c(1, 1, 0, 1), 2, 2,
                 dimnames = list(NULL, c("(Intercept)", "grptrt"))),
      rhs = "~ grp")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb")
    expect_length(st$beta_k, 2L)
    expect_equal(st$beta_k[2], 0)
  })
})
```

- [ ] **Step 2: Run the test — expect failure (functions not found).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the two new `describe` blocks error (`.dd_tmb_build_tmb_data` / `.dd_tmb_default_starts` not found).

- [ ] **Step 3: Implement both helpers in `R/dd-tmb.R`.**

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
#' Build the TMB data list for MixedDiscounting
#'
#' @param prepared Output from [.dd_tmb_prepare_data()].
#' @param design Output from [.dd_tmb_build_design()].
#' @param equation One of "mazur", "exponential".
#' @param family One of "sltb", "gaussian".
#' @return A list whose names match the C++ `DATA_*` macros: `model`, `y`, `x`,
#'   `subject_id` (0-indexed), `X`, `eqn_type`, `family`, `n_obs`, `n_subjects`.
#' @keywords internal
.dd_tmb_build_tmb_data <- function(prepared, design, equation, family) {
  eqn_type <- switch(equation,
    mazur = 0L,
    exponential = 1L,
    stop("Unknown equation: ", equation, call. = FALSE)
  )
  fam_type <- switch(family,
    sltb = 0L,
    gaussian = 1L,
    stop("Unknown family: ", family, call. = FALSE)
  )

  list(
    model = "MixedDiscounting",
    y = as.numeric(prepared$y),
    x = as.numeric(prepared$x),
    subject_id = as.integer(prepared$subject_id),
    X = as.matrix(design$X),
    eqn_type = eqn_type,
    family = fam_type,
    n_obs = as.integer(prepared$n_obs),
    n_subjects = as.integer(prepared$n_subjects)
  )
}


#' Generate default starting values for the TMB discounting model
#'
#' `beta_k` intercept is data-driven: invert the discounting function at the
#' median `y` observed at the minimum delay. `log_sigma_u = log(0.5)`. The
#' generic auxiliary scalar `log_aux` starts at `log(8)` for SLT-beta
#' (population precision phi) and `log(0.1)` for Gaussian (residual SD).
#'
#' @param prepared Output from [.dd_tmb_prepare_data()].
#' @param design Output from [.dd_tmb_build_design()].
#' @param family One of "sltb", "gaussian".
#' @param equation One of "mazur", "exponential" (for the start inversion).
#' @return A parameters list: `beta_k` (length `ncol(X)`), `log_sigma_u`,
#'   `log_aux`, `u` (matrix `n_subjects` x 1).
#' @keywords internal
.dd_tmb_default_starts <- function(prepared, design, family,
                                   equation = "mazur") {
  y <- prepared$y
  x <- prepared$x
  p <- ncol(design$X)
  n_subjects <- prepared$n_subjects

  x_min <- min(x, na.rm = TRUE)
  mu0 <- stats::median(y[x == x_min], na.rm = TRUE)
  # Guard mu0 into the open interval so the inversion is finite.
  mu0 <- min(max(mu0, 1e-3), 1 - 1e-3)

  if (identical(equation, "exponential")) {
    # mu = exp(-k * x_min)  =>  k = -log(mu) / x_min
    k0 <- -log(mu0) / max(x_min, 1e-8)
  } else {
    # mazur: mu = 1 / (1 + k * x_min)  =>  k = (1/mu - 1) / x_min
    k0 <- (1 / mu0 - 1) / max(x_min, 1e-8)
  }
  if (!is.finite(k0) || k0 <= 0) k0 <- 0.01
  log_k0 <- log(k0)

  beta_k <- rep(0, p)
  beta_k[1] <- log_k0

  log_aux <- if (identical(family, "gaussian")) log(0.1) else log(8)

  list(
    beta_k = beta_k,
    log_sigma_u = log(0.5),
    log_aux = log_aux,
    u = matrix(0, nrow = n_subjects, ncol = 1L)
  )
}
```

- [ ] **Step 4: Run the test — expect pass.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: both new blocks pass (`.dd_tmb_build_tmb_data` 2 `it`, `.dd_tmb_default_starts` 3 `it`).

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): add .dd_tmb_build_tmb_data and .dd_tmb_default_starts

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task P.4: `.dd_tmb_run_optimizer()` + `.expand_bounds()` (verbatim port, family-agnostic)

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append `.expand_bounds`, `.dd_tmb_run_optimizer`)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (add `describe` block; pure-R, no compile)

This is a mechanical PORT of beezdemand `.tmb_run_optimizer` and its helper `.expand_bounds`.
The function is already family-agnostic (it only sees `obj$fn`/`obj$gr`), so nothing in the
body changes; we only rename it to the `.dd_` prefix.

- [ ] **Step 1: Write failing test (optimizer on a trivial quadratic obj stub + bounds expansion).**

Append to `tests/testthat/test-fit_dd_tmb.R`:

```r
describe(".expand_bounds()", {
  it("returns the default for every name when bounds is NULL", {
    out <- .expand_bounds(NULL, c("beta_k", "beta_k", "log_aux"), -Inf)
    expect_equal(out, c(-Inf, -Inf, -Inf))
  })
  it("applies a named bound to all matching positions", {
    out <- .expand_bounds(c(beta_k = -2), c("beta_k", "beta_k", "log_aux"), -Inf)
    expect_equal(unname(out), c(-2, -2, -Inf))
  })
  it("warns on an unknown bound name", {
    expect_warning(
      .expand_bounds(c(nope = 1), c("beta_k"), Inf),
      "unknown parameter"
    )
  })
})

describe(".dd_tmb_run_optimizer()", {
  it("minimizes a quadratic via nlminb and normalizes fields", {
    obj <- list(
      par = c(a = 5, b = -5),
      fn = function(p) sum((p - c(1, 2))^2),
      gr = function(p) 2 * (p - c(1, 2))
    )
    ctrl <- list(optimizer = "nlminb", iter_max = 100, eval_max = 200,
                 rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res <- .dd_tmb_run_optimizer(obj, obj$par, ctrl,
                                 user_specified = character(0), verbose = 0)
    expect_equal(unname(res$opt$par), c(1, 2), tolerance = 1e-5)
    expect_equal(res$opt$convergence, 0L)
    expect_true(is.character(res$opt$message))
  })

  it("returns convergence code 99 and finite-Inf objective on optimizer error", {
    obj <- list(par = c(a = 0), fn = function(p) stop("boom"),
                gr = function(p) p)
    ctrl <- list(optimizer = "nlminb", iter_max = 10, eval_max = 20,
                 rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res <- .dd_tmb_run_optimizer(obj, obj$par, ctrl,
                                 user_specified = character(0), verbose = 0)
    expect_equal(res$opt$convergence, 99L)
    expect_equal(res$opt$objective, Inf)
  })
})
```

- [ ] **Step 2: Run the test — expect failure (functions not found).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the two new blocks error (`.expand_bounds` / `.dd_tmb_run_optimizer` not found).

- [ ] **Step 3: PORT `.expand_bounds` and `.dd_tmb_run_optimizer` into `R/dd-tmb.R`.**

Copy **verbatim** from `/Users/brent/Dropbox/GIT/beezdemand/R/tmb-demand.R`:
- `.expand_bounds` — lines **416–430** (copy unchanged; same name).
- `.tmb_run_optimizer` — lines **454–552** — copy the entire body, with the **only**
  change being the function name `.tmb_run_optimizer` -> `.dd_tmb_run_optimizer`. Change
  nothing else (it already references only `obj$fn`/`obj$gr`/`nlminb`/`stats::optim` and is
  family-agnostic).

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
#' Expand partial optimizer bounds to the full parameter vector
#'
#' Ported verbatim from beezdemand `.expand_bounds`.
#' @keywords internal
.expand_bounds <- function(bounds, par_names, default_val) {
  if (is.null(bounds)) return(rep(default_val, length(par_names)))
  result <- rep(default_val, length(par_names))
  names(result) <- par_names
  for (nm in names(bounds)) {
    idx <- which(par_names == nm)
    if (length(idx) > 0) {
      result[idx] <- bounds[nm]
    } else {
      warning("Bounds specified for unknown parameter '", nm, "' (ignored)",
              call. = FALSE)
    }
  }
  result
}


#' Run a single TMB optimization (nlminb or L-BFGS-B), family-agnostic
#'
#' Ported verbatim from beezdemand `.tmb_run_optimizer` (renamed only). Dispatches
#' to `nlminb` or `optim(method = "L-BFGS-B")` and normalizes the return so
#' downstream code sees identical field names regardless of optimizer.
#'
#' @param obj TMB objective object (with `$fn`, `$gr`, `$par`).
#' @param start Named numeric vector of starting values.
#' @param tmb_control Merged control list.
#' @param user_specified Character vector of fields the user set in tmb_control.
#' @param verbose Integer verbosity.
#' @return list(opt = list(par, objective, convergence, message), warnings).
#' @keywords internal
.dd_tmb_run_optimizer <- function(obj, start, tmb_control, user_specified, verbose) {
  optimizer <- tmb_control$optimizer
  iter_max <- tmb_control$iter_max
  eval_max <- tmb_control$eval_max
  rel_tol <- tmb_control$rel_tol

  trace <- if ("trace" %in% user_specified) {
    as.integer(tmb_control$trace)
  } else if (verbose >= 2) {
    1L
  } else {
    0L
  }

  par_names <- names(start)
  lower <- .expand_bounds(tmb_control$lower, par_names, -Inf)
  upper <- .expand_bounds(tmb_control$upper, par_names, Inf)

  opt_warnings <- character(0)

  if (optimizer == "nlminb") {
    opt <- tryCatch(
      withCallingHandlers(
        nlminb(
          start = start,
          objective = obj$fn,
          gradient = obj$gr,
          lower = lower,
          upper = upper,
          control = list(
            eval.max = eval_max,
            iter.max = iter_max,
            rel.tol = rel_tol,
            trace = trace
          )
        ),
        warning = function(w) {
          opt_warnings <<- c(opt_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        list(
          par = start,
          objective = Inf,
          convergence = 99L,
          message = conditionMessage(e)
        )
      }
    )
  } else {
    opt <- tryCatch(
      {
        raw <- withCallingHandlers(
          stats::optim(
            par = start,
            fn = obj$fn,
            gr = obj$gr,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper,
            control = list(
              maxit = iter_max,
              trace = trace
            )
          ),
          warning = function(w) {
            opt_warnings <<- c(opt_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )
        list(
          par = raw$par,
          objective = raw$value,
          convergence = raw$convergence,
          message = raw$message %||% "maximum iterations reached"
        )
      },
      error = function(e) {
        list(
          par = start,
          objective = Inf,
          convergence = 99L,
          message = conditionMessage(e)
        )
      }
    )
  }

  if (is.null(opt$message)) opt$message <- "unknown"

  list(opt = opt, warnings = opt_warnings)
}
```

Note: the `%||%` operator is provided by `R/utils-pipe.R` in the package (already imported);
if not exported, add `#' @importFrom rlang %||%` to the roxygen of `R/dd-tmb.R` or rely on
the package's existing `%||%`.

- [ ] **Step 4: Run the test — expect pass.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: `.expand_bounds` (3 `it`) and `.dd_tmb_run_optimizer` (2 `it`) pass.

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): port .dd_tmb_run_optimizer and .expand_bounds (family-agnostic)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task P.5: `.dd_tmb_multi_start()` — 3 start sets + degenerate-optimum GUARD

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append helper + `phi_min`/guard constants)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (add `describe` block; uses the compiled DLL — `skip_on_cran()` + `skip_if_not_installed("TMB")`, and uses the simulator from `helper-dd-sim.R`)

The multi-start runs 3 starting sets (data-driven, low-k, high-k), keeps the lowest **finite**
nll among **non-degenerate** fits. The GUARD rejects a fit whose `log_aux` floor implies
`phi < phi_min = 0.1` (the φ→0/k→∞ collapse from spec §4.8) **or** whose `beta_k` intercept
is implausibly large in magnitude (`|beta_k[1]| > beta_k_abs_max = 20`, i.e. k outside
`[exp(-20), exp(20)]`). Among the **kept** (non-degenerate) fits we take the lowest nll; if
ALL fits are degenerate we fall back to the lowest-nll fit overall and warn.

- [ ] **Step 1: Write failing test for `.dd_tmb_multi_start` (data-driven sane fit + degenerate rejection).**

Append to `tests/testthat/test-fit_dd_tmb.R`:

```r
describe(".dd_tmb_multi_start()", {
  it("returns a finite-nll fit with a sane intercept on clean data", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(11)
    sim <- .simulate_dd_ip_mixed(n_subjects = 30, log_k_pop = log(0.02),
                                 sigma_u = 0.5, phi = 10, family = "sltb",
                                 equation = "mazur", seed = 11)
    prep <- .dd_tmb_prepare_data(sim, "y", "x", "id")
    design <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "sltb")
    starts <- .dd_tmb_default_starts(prep, design, "sltb", "mazur")
    ctrl <- list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
                 rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res <- .dd_tmb_multi_start(tmb_data, starts, ctrl,
                              user_specified = character(0), verbose = 0)
    expect_true(is.finite(res$opt$objective))
    beta0 <- res$opt$par[names(res$opt$par) == "beta_k"][1]
    expect_lt(abs(beta0), 20)               # not the k->inf collapse
    log_aux <- res$opt$par[["log_aux"]]
    expect_gt(exp(log_aux), 0.1)            # phi above the floor
  })

  it("rejects a degenerate phi->0 candidate in favor of a sane one", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # A single boundary-heavy subject: many exact 0/1 IPs that make the
    # SLT likelihood prefer phi->0, k->inf. With the RE prior + guard the
    # kept fit must stay sane.
    set.seed(70)
    sim <- .simulate_dd_ip_mixed(n_subjects = 25, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 8, family = "sltb",
                                 equation = "mazur", seed = 70)
    # Inject a Jarvis-70-style boundary subject (all 0s and 1s).
    bad <- data.frame(
      id = factor("boundary"),
      x = c(7, 30, 180, 365, 730, 1460, 2920),
      y = c(1, 1, 1, 0, 0, 0, 0)
    )
    sim2 <- rbind(
      data.frame(id = as.character(sim$id), x = sim$x, y = sim$y),
      data.frame(id = as.character(bad$id), x = bad$x, y = bad$y)
    )
    prep <- .dd_tmb_prepare_data(sim2, "y", "x", "id")
    design <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "sltb")
    starts <- .dd_tmb_default_starts(prep, design, "sltb", "mazur")
    ctrl <- list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
                 rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0)
    res <- .dd_tmb_multi_start(tmb_data, starts, ctrl,
                              user_specified = character(0), verbose = 0)
    beta0 <- res$opt$par[names(res$opt$par) == "beta_k"][1]
    log_aux <- res$opt$par[["log_aux"]]
    # population k recovered near truth, NOT 414000-style blowup
    expect_lt(exp(beta0), 1)
    expect_gt(exp(log_aux), 0.1)
  })
})
```

- [ ] **Step 2: Run the test — expect failure (`.dd_tmb_multi_start` not found).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_multi_start()` block errors (function not found). (If the DLL is not built, the two `it` skip — acceptable for this RED step; build the DLL via the C++ phase before the GREEN step.)

- [ ] **Step 3: Implement `.dd_tmb_multi_start` in `R/dd-tmb.R`.**

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
# Guard constants for the degenerate phi->0 / k->inf optimum (spec section 4.8).
.dd_phi_min <- 0.1
.dd_beta_k_abs_max <- 20

#' Multi-start optimization for the TMB discounting model with degenerate-optimum guard
#'
#' Builds 3 starting sets (data-driven, low-k, high-k), runs each through
#' [.dd_tmb_run_optimizer()], and selects the lowest **finite** nll among
#' **non-degenerate** fits. A fit is degenerate when its auxiliary scalar implies
#' `phi < .dd_phi_min` (only meaningful for family == sltb) or `|beta_k[1]|` >
#' `.dd_beta_k_abs_max` (k outside `[exp(-20), exp(20)]`). If every fit is
#' degenerate, the lowest-nll fit is returned with a warning.
#'
#' @param tmb_data TMB data list from [.dd_tmb_build_tmb_data()].
#' @param start_values Default starting list from [.dd_tmb_default_starts()].
#' @param tmb_control Merged control list.
#' @param user_specified Character vector of user-set tmb_control fields.
#' @param verbose Integer verbosity.
#' @return list(obj, opt, nll, start_idx, opt_warnings).
#' @keywords internal
.dd_tmb_multi_start <- function(tmb_data, start_values,
                                tmb_control, user_specified, verbose) {
  is_sltb <- identical(tmb_data$family, 0L)

  # 3 starting sets: data-driven, low-k, high-k.
  start_sets <- vector("list", 3L)
  start_sets[[1]] <- start_values

  s2 <- start_values                          # low-k
  s2$beta_k[1] <- start_values$beta_k[1] - 1.5
  s2$log_sigma_u <- log(0.3)
  start_sets[[2]] <- s2

  s3 <- start_values                          # high-k
  s3$beta_k[1] <- start_values$beta_k[1] + 1.5
  s3$log_sigma_u <- log(0.8)
  start_sets[[3]] <- s3

  .is_degenerate <- function(opt) {
    par <- opt$par
    beta0 <- par[names(par) == "beta_k"][1]
    if (!is.finite(beta0) || abs(beta0) > .dd_beta_k_abs_max) return(TRUE)
    if (is_sltb) {
      la <- par[["log_aux"]]
      if (!is.finite(la) || exp(la) < .dd_phi_min) return(TRUE)
    }
    FALSE
  }

  best_kept_nll <- Inf
  best_kept <- NULL
  best_any_nll <- Inf
  best_any <- NULL

  for (s in seq_along(start_sets)) {
    starts_i <- start_sets[[s]]
    result <- tryCatch({
      obj_i <- TMB::MakeADFun(
        data = tmb_data,
        parameters = starts_i,
        random = "u",
        DLL = "beezdiscounting",
        silent = verbose < 2
      )
      opt_res_i <- .dd_tmb_run_optimizer(
        obj_i, obj_i$par, tmb_control, user_specified, verbose
      )
      list(obj = obj_i, opt = opt_res_i$opt, nll = opt_res_i$opt$objective,
           start_idx = s, opt_warnings = opt_res_i$warnings)
    }, error = function(e) {
      if (verbose >= 2) message(sprintf("  Start set %d failed: %s", s, e$message))
      NULL
    })

    if (is.null(result) || !is.finite(result$nll)) next

    if (result$nll < best_any_nll) {
      best_any_nll <- result$nll
      best_any <- result
    }
    if (!.is_degenerate(result$opt) && result$nll < best_kept_nll) {
      best_kept_nll <- result$nll
      best_kept <- result
    }
  }

  if (is.null(best_any)) {
    stop("All starting value sets failed. ",
         "Check data quality or try different start values.", call. = FALSE)
  }

  if (!is.null(best_kept)) {
    best_result <- best_kept
  } else {
    best_result <- best_any
    warning("All multi-start fits hit the degenerate boundary ",
            "(phi -> 0 / k -> Inf); returning the lowest-nll fit. ",
            "Inspect data for boundary-heavy subjects.", call. = FALSE)
  }

  if (verbose >= 1) {
    message(sprintf("  Multi-start: best NLL = %.2f (start set %d of %d)",
                    best_result$nll, best_result$start_idx, length(start_sets)))
  }

  best_result
}
```

- [ ] **Step 4: Run the test — expect pass (DLL must be built first).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_multi_start()` block passes (2 `it` green); the sane-intercept and degenerate-rejection assertions hold.

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): add .dd_tmb_multi_start with degenerate-optimum guard

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task P.6: `.dd_tmb_extract_estimates()` — sdreport + Hessian-PD gate + log_aux rename

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append helper)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (add `describe` block; compiled — `skip_on_cran()` + `skip_if_not_installed("TMB")`)

This is a SIMPLIFIED port of beezdemand `.tmb_extract_estimates` (lines **722–824**). Single
random column per subject (`re_dim_total = 1`); no vectorized `logsigma`/`rho_raw`. The
auxiliary scalar `log_aux` is **renamed** in the returned `coefficients`/`se` to `log_phi`
(sltb) or `log_sigma_e` (gaussian).

- [ ] **Step 1: Write failing test for `.dd_tmb_extract_estimates`.**

Append to `tests/testthat/test-fit_dd_tmb.R`:

```r
describe(".dd_tmb_extract_estimates()", {
  it("returns sane coefficients, renames log_aux, and gates pdHess", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(21)
    sim <- .simulate_dd_ip_mixed(n_subjects = 30, log_k_pop = log(0.02),
                                 sigma_u = 0.5, phi = 10, family = "sltb",
                                 equation = "mazur", seed = 21)
    prep <- .dd_tmb_prepare_data(sim, "y", "x", "id")
    design <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "sltb")
    starts <- .dd_tmb_default_starts(prep, design, "sltb", "mazur")
    obj <- TMB::MakeADFun(tmb_data, starts, random = "u",
                          DLL = "beezdiscounting", silent = TRUE)
    opt_res <- .dd_tmb_run_optimizer(
      obj, obj$par,
      list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
           rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0),
      character(0), 0)
    est <- .dd_tmb_extract_estimates(obj, opt_res$opt,
                                     n_subjects = prep$n_subjects,
                                     family = "sltb", verbose = 0)
    expect_true("log_phi" %in% names(est$coefficients))
    expect_false("log_aux" %in% names(est$coefficients))
    expect_true("log_phi" %in% names(est$se))
    expect_equal(dim(est$u_hat), c(prep$n_subjects, 1L))
    expect_true(is.logical(est$hessian_pd))
  })

  it("renames log_aux to log_sigma_e under gaussian", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(22)
    sim <- .simulate_dd_ip_mixed(n_subjects = 25, log_k_pop = log(0.02),
                                 sigma_u = 0.5, sigma_e = 0.08,
                                 family = "gaussian", equation = "mazur",
                                 seed = 22)
    prep <- .dd_tmb_prepare_data(sim, "y", "x", "id")
    design <- .dd_tmb_build_design(prep$data)
    tmb_data <- .dd_tmb_build_tmb_data(prep, design, "mazur", "gaussian")
    starts <- .dd_tmb_default_starts(prep, design, "gaussian", "mazur")
    obj <- TMB::MakeADFun(tmb_data, starts, random = "u",
                          DLL = "beezdiscounting", silent = TRUE)
    opt_res <- .dd_tmb_run_optimizer(
      obj, obj$par,
      list(optimizer = "nlminb", iter_max = 1000, eval_max = 2000,
           rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0),
      character(0), 0)
    est <- .dd_tmb_extract_estimates(obj, opt_res$opt,
                                     n_subjects = prep$n_subjects,
                                     family = "gaussian", verbose = 0)
    expect_true("log_sigma_e" %in% names(est$coefficients))
  })
})
```

- [ ] **Step 2: Run the test — expect failure (function not found / skipped if no DLL).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_extract_estimates()` block errors (function not found).

- [ ] **Step 3: Implement `.dd_tmb_extract_estimates` in `R/dd-tmb.R`.**

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
#' Extract estimates from a TMB discounting fit
#'
#' Computes `TMB::sdreport`, gates on `isTRUE(sdr$pdHess)` (warn, never abort),
#' and renames the generic auxiliary scalar `log_aux` to `log_phi` (sltb) or
#' `log_sigma_e` (gaussian) in the returned `coefficients`/`se`.
#'
#' @param obj TMB objective object.
#' @param opt Normalized optimizer result (from [.dd_tmb_run_optimizer()]).
#' @param n_subjects Integer number of subjects.
#' @param family One of "sltb", "gaussian".
#' @param verbose Integer verbosity.
#' @return list(coefficients, se, sdr, variance_components, u_hat, hessian_pd).
#' @keywords internal
.dd_tmb_extract_estimates <- function(obj, opt, n_subjects, family, verbose = 1) {
  sdr <- tryCatch(
    TMB::sdreport(obj),
    error = function(e1) {
      sdr2 <- tryCatch(
        TMB::sdreport(obj, getJointPrecision = FALSE),
        error = function(e2) NULL
      )
      if (is.null(sdr2) && verbose >= 1) {
        warning("Standard error computation failed: ", e1$message)
      }
      sdr2
    }
  )

  hessian_pd <- NA
  if (!is.null(sdr)) {
    hessian_pd <- isTRUE(sdr$pdHess)
    if (!hessian_pd && verbose >= 1) {
      cli::cli_warn(c(
        "!" = "Hessian is not positive definite ({.code pdHess = FALSE}).",
        "i" = "Standard errors, p-values, and confidence intervals may be unreliable.",
        "i" = "Consider simplifying the model or checking data quality."
      ))
    }
  }

  par_full <- opt$par
  par_names <- names(par_full)

  coefficients <- par_full
  se_vec <- rep(NA_real_, length(par_full))
  names(se_vec) <- par_names

  if (!is.null(sdr)) {
    fixed_summary <- summary(sdr, "fixed")
    .fill_vector_se <- function(name) {
      idx <- which(par_names == name)
      if (length(idx) == 0L) return(invisible(NULL))
      rows <- fixed_summary[rownames(fixed_summary) == name, , drop = FALSE]
      if (nrow(rows) == length(idx)) se_vec[idx] <<- rows[, "Std. Error"]
    }
    .fill_vector_se("beta_k")
    .fill_vector_se("log_sigma_u")
    .fill_vector_se("log_aux")

    re_summary <- tryCatch(summary(sdr, "random"), error = function(e) NULL)
    if (!is.null(re_summary)) {
      u_hat <- matrix(re_summary[, "Estimate"], nrow = n_subjects, ncol = 1L)
    } else {
      u_hat <- matrix(0, nrow = n_subjects, ncol = 1L)
    }
  } else {
    u_hat <- matrix(0, nrow = n_subjects, ncol = 1L)
  }

  variance_components <- NULL
  if (!is.null(sdr)) {
    variance_components <- tryCatch(summary(sdr, "report"), error = function(e) NULL)
  }

  # Rename the generic auxiliary scalar.
  aux_name <- if (identical(family, "gaussian")) "log_sigma_e" else "log_phi"
  names(coefficients)[names(coefficients) == "log_aux"] <- aux_name
  names(se_vec)[names(se_vec) == "log_aux"] <- aux_name

  list(
    coefficients = coefficients,
    se = se_vec,
    sdr = sdr,
    variance_components = variance_components,
    u_hat = u_hat,
    hessian_pd = hessian_pd
  )
}
```

- [ ] **Step 4: Run the test — expect pass.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_extract_estimates()` block passes (2 `it` green).

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): add .dd_tmb_extract_estimates (sdreport, pdHess gate, log_aux rename)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task P.7: `.dd_tmb_compute_subject_pars()` — `k_i = exp(Xbeta + sigma_u * u_i)`

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append helper)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (add `describe` block; pure-R, no compile)

Single random intercept: `re_i = sigma_u * u_i` (non-centered, matching the C++ predictor
`log_k_i = X.row(i)*beta_k + sigma_u * u(subj,0)`). For the MVP intercept-only design,
`Xbeta` is the same for all subjects (the intercept); subject k differs only via `u_i`.
Output `data.frame(id, u_i, k [, phi])` — include a `phi` column only for sltb.

- [ ] **Step 1: Write failing test for `.dd_tmb_compute_subject_pars`.**

Append to `tests/testthat/test-fit_dd_tmb.R`:

```r
describe(".dd_tmb_compute_subject_pars()", {
  it("computes k_i = exp(beta0 + sigma_u * u_i) for an intercept-only fit (sltb)", {
    coefs <- c(beta_k = log(0.02), log_sigma_u = log(0.5), log_phi = log(8))
    u_hat <- matrix(c(-1, 0, 2), ncol = 1L)
    sp <- .dd_tmb_compute_subject_pars(
      coefficients = coefs, u_hat = u_hat,
      subject_levels = c("a", "b", "c"),
      equation = "mazur", family = "sltb")
    expect_named(sp, c("id", "u_i", "k", "phi"))
    sigma_u <- 0.5
    expect_equal(sp$k, exp(log(0.02) + sigma_u * c(-1, 0, 2)), tolerance = 1e-10)
    expect_equal(sp$u_i, c(-1, 0, 2))
    expect_equal(unique(sp$phi), 8, tolerance = 1e-10)
  })

  it("omits the phi column for gaussian fits", {
    coefs <- c(beta_k = log(0.02), log_sigma_u = log(0.5), log_sigma_e = log(0.1))
    u_hat <- matrix(c(0, 1), ncol = 1L)
    sp <- .dd_tmb_compute_subject_pars(
      coefficients = coefs, u_hat = u_hat,
      subject_levels = c("a", "b"),
      equation = "mazur", family = "gaussian")
    expect_named(sp, c("id", "u_i", "k"))
    expect_false("phi" %in% names(sp))
  })
})
```

- [ ] **Step 2: Run the test — expect failure (function not found).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_compute_subject_pars()` block errors (function not found).

- [ ] **Step 3: Implement `.dd_tmb_compute_subject_pars` in `R/dd-tmb.R`.**

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
#' Compute subject-specific discounting parameters
#'
#' Reconstructs each subject's `k_i = exp(beta_k[1] + sigma_u * u_i)` using the
#' non-centered predictor that matches the C++ template
#' (`log_k_i = X.row(i)*beta_k + sigma_u * u(subj,0)`). For the intercept-only
#' MVP, `Xbeta` equals the intercept for all subjects; subject k differs via
#' `u_i`. A `phi` column (population precision) is added only for `family =
#' "sltb"`.
#'
#' @param coefficients Named coefficient vector (with `beta_k`, `log_sigma_u`,
#'   and `log_phi` or `log_sigma_e`).
#' @param u_hat Matrix `n_subjects` x 1 of standardized random effects.
#' @param subject_levels Character vector of subject ids (length n_subjects).
#' @param equation One of "mazur", "exponential" (reserved; k is equation-free).
#' @param family One of "sltb", "gaussian".
#' @return data.frame(id, u_i, k[, phi]).
#' @keywords internal
.dd_tmb_compute_subject_pars <- function(coefficients, u_hat, subject_levels,
                                         equation, family) {
  beta_k <- unname(coefficients[names(coefficients) == "beta_k"])
  beta0 <- beta_k[1]
  sigma_u <- exp(unname(coefficients[["log_sigma_u"]]))
  u_i <- as.numeric(u_hat[, 1L])

  log_k_i <- beta0 + sigma_u * u_i
  k_i <- exp(log_k_i)

  out <- data.frame(
    id = subject_levels,
    u_i = u_i,
    k = k_i,
    stringsAsFactors = FALSE
  )

  if (identical(family, "sltb")) {
    out$phi <- exp(unname(coefficients[["log_phi"]]))
  }

  out
}
```

- [ ] **Step 4: Run the test — expect pass.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `.dd_tmb_compute_subject_pars()` block passes (2 `it` green).

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): add .dd_tmb_compute_subject_pars (non-centered k_i)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task P.8: `fit_dd_tmb()` — full pipeline wiring + recovery & degenerate regression tests

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R` (append the public `fit_dd_tmb`)
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/NAMESPACE` (via roxygen `@export`)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd_tmb.R` (add the public-API `describe` block)

`fit_dd_tmb()` calls `.dd_validate_ip` -> `.dd_normalize_re` -> `.dd_tmb_prepare_data` ->
`.dd_tmb_build_design` -> `.dd_tmb_build_tmb_data` -> `.dd_tmb_default_starts` -> (multi-start
or single `TMB::MakeADFun(DLL="beezdiscounting", random="u")` + `.dd_tmb_run_optimizer`) ->
`.dd_tmb_extract_estimates` -> `.dd_tmb_compute_subject_pars`, then assembles the
`beezdiscounting_tmb` object per the pinned contract.

- [ ] **Step 1: Write failing recovery + degenerate-regression tests for `fit_dd_tmb`.**

Append to `tests/testthat/test-fit_dd_tmb.R`:

```r
describe("fit_dd_tmb() recovery", {
  it("recovers population k within 0.15 (sltb x mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(101)
    sim <- .simulate_dd_ip_mixed(n_subjects = 60, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 10, family = "sltb",
                                 equation = "mazur", seed = 101)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb", verbose = 0)
    expect_s3_class(fit, "beezdiscounting_tmb")
    beta0 <- fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1]
    expect_equal(exp(beta0), 0.01, tolerance = 0.15)
    expect_true(fit$converged)
    expect_true(fit$se_available)
  })

  it("recovers population k within 0.15 (gaussian x mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(102)
    sim <- .simulate_dd_ip_mixed(n_subjects = 60, log_k_pop = log(0.01),
                                 sigma_u = 0.6, sigma_e = 0.08,
                                 family = "gaussian", equation = "mazur",
                                 seed = 102)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "gaussian", verbose = 0)
    beta0 <- fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1]
    expect_equal(exp(beta0), 0.01, tolerance = 0.15)
    expect_true(fit$converged)
  })

  it("recovers k for the exponential equation (sltb, looser tol)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(103)
    sim <- .simulate_dd_ip_mixed(n_subjects = 60, log_k_pop = log(0.005),
                                 sigma_u = 0.5, phi = 12, family = "sltb",
                                 equation = "exponential", seed = 103)
    fit <- fit_dd_tmb(sim, equation = "exponential", family = "sltb",
                      verbose = 0)
    expect_true(fit$converged)
    # correlation of subject log-k recovery over the grid (looser per contract)
    truth_k <- 0.005
    beta0 <- fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1]
    expect_equal(exp(beta0), truth_k, tolerance = 0.30)
  })

  it("recovers k for gaussian x exponential", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(104)
    sim <- .simulate_dd_ip_mixed(n_subjects = 60, log_k_pop = log(0.005),
                                 sigma_u = 0.5, sigma_e = 0.06,
                                 family = "gaussian", equation = "exponential",
                                 seed = 104)
    fit <- fit_dd_tmb(sim, equation = "exponential", family = "gaussian",
                      verbose = 0)
    expect_true(fit$converged)
    beta0 <- fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1]
    expect_equal(exp(beta0), 0.005, tolerance = 0.30)
  })
})

describe("fit_dd_tmb() object shape", {
  it("assembles the contract fields and subject_pars", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(105)
    sim <- .simulate_dd_ip_mixed(n_subjects = 40, family = "sltb",
                                 equation = "mazur", seed = 105)
    fit <- fit_dd_tmb(sim, verbose = 0)
    expect_s3_class(fit, "beezdiscounting_tmb")
    expect_true(all(c("opt", "model", "sdr", "param_info", "formula_details",
                      "subject_pars", "loglik", "AIC", "BIC", "converged",
                      "se_available", "data", "data_all",
                      "coercion_info") %in% names(fit)))
    expect_named(fit$subject_pars, c("id", "u_i", "k", "phi"))
    expect_equal(nrow(fit$subject_pars), 40L)
    expect_equal(fit$param_info$family, "sltb")
    expect_equal(fit$param_info$equation, "mazur")
    expect_equal(fit$param_info$n_random_effects, 1L)
    expect_equal(fit$param_info$id_var, "id")
    expect_true("log_phi" %in% names(fit$model$coefficients))
  })

  it("keeps a sane population k on a boundary-heavy dataset (regression)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    set.seed(70)
    sim <- .simulate_dd_ip_mixed(n_subjects = 25, log_k_pop = log(0.01),
                                 sigma_u = 0.6, phi = 8, family = "sltb",
                                 equation = "mazur", seed = 70)
    bad <- data.frame(
      id = "boundary",
      x = c(7, 30, 180, 365, 730, 1460, 2920),
      y = c(1, 1, 1, 0, 0, 0, 0)
    )
    sim2 <- rbind(
      data.frame(id = as.character(sim$id), x = sim$x, y = sim$y), bad)
    fit <- fit_dd_tmb(sim2, equation = "mazur", family = "sltb", verbose = 0)
    beta0 <- fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1]
    # NOT the k -> 414000 collapse; population k stays plausible
    expect_lt(exp(beta0), 1)
    expect_gt(exp(fit$model$coefficients[["log_phi"]]), 0.1)
  })
})
```

- [ ] **Step 2: Run the tests — expect failure (`fit_dd_tmb` not found / skipped if no DLL).**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: the `fit_dd_tmb()` blocks error (function not found).

- [ ] **Step 3: Implement `fit_dd_tmb()` in `R/dd-tmb.R`.**

Append to `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb.R`:

```r
#' Fit an indifference-point mixed-effects discounting model via TMB
#'
#' Fits a 1-parameter discounting model (Mazur hyperbolic or exponential) with a
#' random intercept on `log k`, between-subject fixed effects, and either an
#' SLT-beta or Gaussian observation family, using Template Model Builder for
#' exact AD + Laplace approximation.
#'
#' @param data Long data frame with subject id, delay, and indifference
#'   proportion columns.
#' @param y_var,x_var,id_var Column names (defaults `"y"`, `"x"`, `"id"`).
#' @param equation One of `"mazur"`, `"exponential"`.
#' @param family Observation family: `"sltb"` (default) or `"gaussian"`.
#' @param random_effects RE formula (MVP: `k ~ 1`).
#' @param factors Character vector of between-subject factor names.
#' @param factor_interaction Logical; include a pairwise factor interaction.
#' @param continuous_covariates Character vector of covariate names.
#' @param ll Optional larger-later reward for `amount`-scale coercion.
#' @param response_scale One of `"proportion"`, `"percent"`, `"amount"`.
#' @param start_values Optional named list overriding defaults.
#' @param tmb_control Optimizer control list.
#' @param multi_start Logical; if `TRUE` (default), run the 3-set guarded
#'   multi-start.
#' @param verbose Integer verbosity (0 silent, 1 progress, 2 debug).
#' @param ... Reserved.
#' @return An object of class `beezdiscounting_tmb`.
#' @export
fit_dd_tmb <- function(data,
                       y_var = "y", x_var = "x", id_var = "id",
                       equation = c("mazur", "exponential"),
                       family = c("sltb", "gaussian"),
                       random_effects = k ~ 1,
                       factors = NULL,
                       factor_interaction = FALSE,
                       continuous_covariates = NULL,
                       ll = NULL,
                       response_scale = c("proportion", "percent", "amount"),
                       start_values = NULL,
                       tmb_control = list(iter_max = 1000, eval_max = 2000),
                       multi_start = TRUE,
                       verbose = 1,
                       ...) {
  cl <- match.call()
  equation <- match.arg(equation)
  family <- match.arg(family)
  response_scale <- match.arg(response_scale)

  # 1. Validate + coerce/clamp (warns; errors on missing/NA y).
  validated <- .dd_validate_ip(data, y_var = y_var, x_var = x_var,
                               id_var = id_var, ll = ll,
                               response_scale = response_scale)
  long <- validated$data            # canonical id/x/y
  coercion_info <- validated$coercion_info

  # 2. Random-effects normalization (MVP: single intercept-only block).
  re_norm <- .dd_normalize_re(random_effects, data = long)
  n_random_effects <- 1L

  # 3. Carry factor / covariate columns onto the cleaned long frame so the
  #    design matrix can be built on the same rows the model is fit on.
  extra_cols <- unique(c(factors, continuous_covariates))
  extra_cols <- intersect(extra_cols, names(data))
  if (length(extra_cols) > 0L) {
    # data and `long` are row-aligned only when no rows were dropped; rejoin
    # by the validator's preserved order via the original row index.
    long <- cbind(long, data[seq_len(nrow(long)), extra_cols, drop = FALSE])
    names(long)[(ncol(long) - length(extra_cols) + 1L):ncol(long)] <- extra_cols
  }

  # 4. Prepare data (0-indexed subject_id, NA drop, row coherence).
  prepared <- .dd_tmb_prepare_data(long, y_var = "y", x_var = "x",
                                   id_var = "id")

  # Restrict the design frame to the surviving rows (prepared$data is cleaned).
  design_data <- prepared$data
  if (length(extra_cols) > 0L) {
    design_data <- cbind(design_data,
                         long[match(rownames(prepared$data), rownames(long)),
                              extra_cols, drop = FALSE])
    # match() may fail if rownames differ; fall back to direct attach when
    # the validator preserved order and no NA rows were dropped.
    if (anyNA(design_data[, extra_cols, drop = FALSE]) &&
        nrow(long) == nrow(prepared$data)) {
      design_data <- cbind(prepared$data, long[, extra_cols, drop = FALSE])
    }
  }

  # 5. Fixed-effect design for log k.
  design <- .dd_tmb_build_design(design_data, factors = factors,
                                 factor_interaction = factor_interaction,
                                 continuous_covariates = continuous_covariates)

  # 6. TMB data + default starts.
  tmb_data <- .dd_tmb_build_tmb_data(prepared, design, equation, family)
  default_starts <- .dd_tmb_default_starts(prepared, design, family, equation)
  if (!is.null(start_values)) {
    for (nm in names(start_values)) {
      if (nm %in% names(default_starts)) default_starts[[nm]] <- start_values[[nm]]
    }
  }

  # 7. Merge control defaults.
  default_control <- list(
    iter_max = 1000, eval_max = 2000, optimizer = "nlminb",
    rel_tol = 1e-10, lower = NULL, upper = NULL, trace = 0
  )
  user_specified <- names(tmb_control)
  tmb_control <- utils::modifyList(default_control, tmb_control)

  if (verbose >= 1) {
    message(sprintf("Fitting TMB mixed-effects discounting model (%s, %s)...",
                    equation, family))
    message(sprintf("  Subjects: %d, Observations: %d",
                    prepared$n_subjects, prepared$n_obs))
  }

  # 8. Optimize (multi-start or single).
  opt_warnings <- character(0)
  if (isTRUE(multi_start)) {
    result <- .dd_tmb_multi_start(tmb_data, default_starts, tmb_control,
                                  user_specified, verbose)
    obj <- result$obj
    opt <- result$opt
    opt_warnings <- result$opt_warnings %||% character(0)
  } else {
    obj <- TMB::MakeADFun(tmb_data, default_starts, random = "u",
                          DLL = "beezdiscounting", silent = verbose < 2)
    opt_res <- .dd_tmb_run_optimizer(obj, obj$par, tmb_control,
                                     user_specified, verbose)
    opt <- opt_res$opt
    opt_warnings <- opt_res$warnings
  }

  converged <- isTRUE(opt$convergence == 0)
  try(obj$fn(opt$par), silent = TRUE)

  # 9. Extract estimates (sdreport, pdHess gate, log_aux rename).
  estimates <- .dd_tmb_extract_estimates(obj, opt,
                                         n_subjects = prepared$n_subjects,
                                         family = family, verbose = verbose)

  # 10. Subject-specific parameters.
  subject_pars <- .dd_tmb_compute_subject_pars(
    coefficients = estimates$coefficients,
    u_hat = estimates$u_hat,
    subject_levels = prepared$subject_levels,
    equation = equation,
    family = family
  )

  # 11. Likelihood / IC.
  nll <- opt$objective
  loglik <- -nll
  n_fixed_params <- length(opt$par)
  aic <- 2 * nll + 2 * n_fixed_params
  bic <- 2 * nll + n_fixed_params * log(prepared$n_obs)

  has_phi <- identical(family, "sltb")

  result_obj <- structure(
    list(
      call = cl,
      opt = opt,
      model = list(
        coefficients = estimates$coefficients,
        se = estimates$se,
        variance_components = estimates$variance_components
      ),
      sdr = estimates$sdr,
      hessian_pd = estimates$hessian_pd,
      param_info = list(
        equation = equation,
        family = family,
        has_phi = has_phi,
        n_obs = prepared$n_obs,
        n_subjects = prepared$n_subjects,
        n_random_effects = n_random_effects,
        subject_levels = prepared$subject_levels,
        id_var = id_var,
        x_var = x_var,
        y_var = y_var,
        random_effects_parsed = re_norm
      ),
      formula_details = list(X = design$X, rhs = design$rhs),
      subject_pars = subject_pars,
      loglik = loglik,
      AIC = aic,
      BIC = bic,
      converged = converged,
      se_available = !is.null(estimates$sdr),
      opt_warnings = opt_warnings,
      data = prepared$data,
      data_all = long,
      coercion_info = coercion_info
    ),
    class = "beezdiscounting_tmb"
  )

  if (verbose >= 1) {
    if (converged) {
      message(sprintf("  Converged (NLL = %.2f). Done.", nll))
    } else {
      message(sprintf("  WARNING: Did not converge (code %s: %s).",
                      opt$convergence, opt$message))
    }
  }

  result_obj
}
```

- [ ] **Step 4: Document (regenerate NAMESPACE) and run the tests — expect pass.**

Run command:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::document(quiet=TRUE); devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd_tmb.R", reporter="summary")'
```
Expected output: NAMESPACE now contains `export(fit_dd_tmb)`; all `describe` blocks in
`test-fit_dd_tmb.R` pass — recovery within 0.15 (mazur) / 0.30 (exponential) for both
families, `converged`/`se_available` TRUE, the object-shape assertions hold, and the
boundary-heavy regression keeps `exp(beta0) < 1` and `phi > 0.1`.

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && git add R/dd-tmb.R NAMESPACE man/fit_dd_tmb.Rd tests/testthat/test-fit_dd_tmb.R && git commit -m "feat(tmb): add public fit_dd_tmb() wiring the dd-tmb pipeline

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
## Phase M — S3 methods

Implements the broom/stats S3 surface for class `beezdiscounting_tmb` in
`R/dd-tmb-methods.R`, mirroring the contracts in
`/Users/brent/Dropbox/GIT/beezdemand/R/tmb-methods.R`. Every method consumes the
fit-object structure from the PINNED CONTRACT: `object$model$coefficients`
(named numeric incl. `beta_k`, `log_sigma_u`, and the renamed `log_phi` /
`log_sigma_e`), `object$model$se` (parallel named SEs — added here as a thin
accessor over the sdreport), `object$sdr`, `object$subject_pars`
(data.frame `id`, `u_i`, `k` [, `phi`]), `object$param_info`,
`object$formula_details$X`, `object$loglik`, `object$AIC`, `object$BIC`,
`object$converged`, `object$se_available`, and `object$data`.

This phase assumes Phase P (param-space) has shipped `.dd_transform_coef_table()`
and `.dd_transform_est_se()` (ported from beezdemand, regex retargeted to
`^k($|_)|^s($|_)|^phi($|_)`), and that Phase F (`fit_dd_tmb`) produces the fit
object. Where a helper is shared across methods (the term-name builder, the
fitted/resid back end), it is defined once here.

---

### Task M.0: Shared term-name + SE helpers

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing test for `.dd_tmb_build_term_names()` and `.dd_tmb_model_se()`**

Append to `tests/testthat/test-dd-tmb-methods.R` (create file with the header
below first):

```r
# tests/testthat/test-dd-tmb-methods.R
skip_on_cran()
skip_if_not_installed("TMB")

# A small SLT-beta Mazur fit reused by every method test in this file.
fit_for_methods <- function() {
  set.seed(101)
  dat <- .simulate_dd_ip_mixed(
    n_subjects = 30, family = "sltb", equation = "mazur",
    log_k_pop = log(0.01), sigma_u = 0.6, phi = 10, seed = 101
  )
  fit_dd_tmb(dat, equation = "mazur", family = "sltb",
             random_effects = k ~ 1, multi_start = TRUE, verbose = 0)
}

describe(".dd_tmb_build_term_names()", {
  it("maps beta_k columns to k:<design colname> and passes aux through", {
    fit <- fit_for_methods()
    tn <- .dd_tmb_build_term_names(fit)
    expect_true("k:(Intercept)" %in% tn$term)
    # aux + sigma_u rows keep their raw display names
    expect_true(any(tn$term %in% c("log_phi", "log_sigma_e")))
    expect_true("log_sigma_u" %in% tn$term)
    expect_length(tn$k_idx, ncol(fit$formula_details$X))
    expect_type(tn$other_idx, "integer")
  })
})

describe(".dd_tmb_model_se()", {
  it("returns a named SE vector aligned to coefficients", {
    fit <- fit_for_methods()
    se <- .dd_tmb_model_se(fit)
    co <- fit$model$coefficients
    expect_named(se, names(co))
    expect_length(se, length(co))
  })
})
```

- [ ] **Step 2: run (expect fail)**

Run command:
```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected output: errors — `could not find function ".dd_tmb_build_term_names"`.

- [ ] **Step 3: implement the file header + helpers**

Create `R/dd-tmb-methods.R`:

```r
# ==============================================================================
# S3 Methods for beezdiscounting_tmb Objects
# ==============================================================================

#' Build display term names from a beezdiscounting_tmb coefficient vector
#'
#' Maps raw optimizer names (`beta_k`, `log_sigma_u`, `log_phi` / `log_sigma_e`)
#' to readable display names. `beta_k` columns become `k:<design colname>`
#' using `formula_details$X` colnames; every other coefficient keeps its raw
#' name. Shared by `tidy()`, `summary()`, `confint()`.
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param nms Character vector of raw parameter names (default
#'   `names(object$model$coefficients)`).
#' @return List with `term` (display names), `k_idx` (beta_k positions),
#'   `other_idx` (non-beta positions).
#' @keywords internal
.dd_tmb_build_term_names <- function(object, nms = NULL) {
  if (is.null(nms)) {
    nms <- names(object$model$coefficients)
  }
  k_idx <- which(nms == "beta_k")
  other_idx <- which(nms != "beta_k")

  k_colnames <- colnames(object$formula_details$X)
  if (is.null(k_colnames)) {
    k_colnames <- paste0("X", seq_along(k_idx))
  }

  term <- character(length(nms))
  term[k_idx] <- paste0("k:", k_colnames)
  term[other_idx] <- nms[other_idx]

  list(
    term = term,
    k_idx = k_idx,
    other_idx = as.integer(other_idx)
  )
}


#' Named standard-error vector aligned to the coefficient vector
#'
#' Pulls fixed-effect SEs from the sdreport. The optimizer parameter vector
#' (`beta_k`, `log_sigma_u`, `log_aux`) is what `sdreport$par.fixed` /
#' `sdreport$cov.fixed` cover; entries are renamed to match
#' `names(object$model$coefficients)` (i.e. `log_aux` already displayed as
#' `log_phi` / `log_sigma_e` in the fit object). Returns `NA` SEs when the
#' sdreport is unavailable.
#'
#' @param object A `beezdiscounting_tmb` object.
#' @return Named numeric vector parallel to `object$model$coefficients`.
#' @keywords internal
.dd_tmb_model_se <- function(object) {
  co <- object$model$coefficients
  if (!is.null(object$model$se)) {
    se <- object$model$se
    # Defensive: align to coefficient order/names.
    return(stats::setNames(unname(se)[match(names(co), names(se))], names(co)))
  }
  sdr <- object$sdr
  if (is.null(sdr) || isFALSE(object$se_available)) {
    return(stats::setNames(rep(NA_real_, length(co)), names(co)))
  }
  sd_fixed <- sqrt(diag(as.matrix(sdr$cov.fixed)))
  raw_nms <- names(sdr$par.fixed)
  # `log_aux` is the optimizer name; the fit object renames it to log_phi /
  # log_sigma_e. Translate so the SE vector lines up with coefficients.
  aux_name <- intersect(c("log_phi", "log_sigma_e"), names(co))
  raw_nms[raw_nms == "log_aux"] <- if (length(aux_name)) aux_name[1] else "log_aux"
  se <- stats::setNames(rep(NA_real_, length(co)), names(co))
  hit <- match(names(co), raw_nms)
  se[!is.na(hit)] <- sd_fixed[hit[!is.na(hit)]]
  se
}
```

- [ ] **Step 4: run (expect pass)**

Run command:
```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected output: 2 `describe` blocks pass; 0 failures.

- [ ] **Step 5: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R tests/testthat/test-dd-tmb-methods.R
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): term-name + SE helpers for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.1: `logLik` / `AIC` / `BIC` / `nobs`

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("logLik / AIC / BIC / nobs", {
  fit <- fit_for_methods()

  it("logLik carries df and nobs attributes", {
    ll <- logLik(fit)
    expect_s3_class(ll, "logLik")
    expect_equal(as.numeric(ll), fit$loglik)
    expect_equal(attr(ll, "df"), length(fit$opt$par))
    expect_equal(attr(ll, "nobs"), fit$param_info$n_obs)
  })

  it("AIC matches the stored value and honours custom k", {
    expect_equal(AIC(fit), fit$AIC)
    np <- length(fit$opt$par)
    expect_equal(AIC(fit, k = 3), 2 * (-fit$loglik) + 3 * np)
  })

  it("BIC and nobs match stored values", {
    expect_equal(BIC(fit), fit$BIC)
    expect_equal(nobs(fit), fit$param_info$n_obs)
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — `no applicable method for 'logLik' applied to an object of class "beezdiscounting_tmb"`.

- [ ] **Step 3: implement**

Append to `R/dd-tmb-methods.R`:

```r
# --- logLik / AIC / BIC / nobs ---

#' @export
logLik.beezdiscounting_tmb <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- length(object$opt$par)
  attr(ll, "nobs") <- object$param_info$n_obs
  class(ll) <- "logLik"
  ll
}

#' @export
AIC.beezdiscounting_tmb <- function(object, ..., k = 2) {
  if (k != 2) {
    nll <- -object$loglik
    n_params <- length(object$opt$par)
    return(2 * nll + k * n_params)
  }
  object$AIC
}

#' @export
BIC.beezdiscounting_tmb <- function(object, ...) {
  object$BIC
}

#' @export
nobs.beezdiscounting_tmb <- function(object, ...) {
  object$param_info$n_obs %||% nrow(object$data)
}
```

(`%||%` is exported from `R/utils.R` / `R/utils-pipe.R` already in the package;
no new import needed.)

- [ ] **Step 4: add roxygen `@export`s to NAMESPACE**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting")'
```
Expected: NAMESPACE gains `S3method(AIC,beezdiscounting_tmb)`,
`S3method(BIC,beezdiscounting_tmb)`, `S3method(logLik,beezdiscounting_tmb)`,
`S3method(nobs,beezdiscounting_tmb)`.

- [ ] **Step 5: run (expect pass)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: `logLik / AIC / BIC / nobs` block passes; 0 failures.

- [ ] **Step 6: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): logLik/AIC/BIC/nobs for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.2: `coef` / `fixef` / `ranef`

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("coef / fixef", {
  fit <- fit_for_methods()

  it("coef() returns the raw named optimizer vector", {
    co <- coef(fit)
    expect_type(co, "double")
    expect_named(co)
    expect_true("beta_k" %in% names(co))
    expect_true("log_sigma_u" %in% names(co))
    expect_true(any(c("log_phi", "log_sigma_e") %in% names(co)))
    expect_identical(co, fit$model$coefficients)
  })

  it("fixef() returns the same named numeric vector as coef()", {
    expect_identical(nlme::fixef(fit), coef(fit))
  })
})

describe("ranef", {
  fit <- fit_for_methods()

  it("returns id, u_i, and per-subject k for every subject", {
    re <- ranef(fit)
    expect_s3_class(re, "data.frame")
    expect_named(re, c("id", "u_i", "k"))
    expect_equal(nrow(re), fit$param_info$n_subjects)
    expect_equal(re$id, fit$subject_pars$id)
    expect_equal(re$u_i, fit$subject_pars$u_i)
    expect_true(all(re$k > 0))
  })

  it("includes phi when family = sltb has subject phi column", {
    # MVP: phi is population-level, so subject_pars has no phi column and
    # ranef() exposes only u_i and k. Guard the contract explicitly.
    re <- ranef(fit)
    expect_false("phi" %in% names(re))
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — no `coef` / `fixef` / `ranef` method for class.

- [ ] **Step 3: implement**

Append to `R/dd-tmb-methods.R`:

```r
# --- coef / fixef / ranef ---

#' Extract coefficients from a TMB discounting model
#'
#' Returns the optimizer's flat named parameter vector: `beta_k` (one entry per
#' fixed-effect design column, on the log-k scale), `log_sigma_u`, and the
#' auxiliary parameter (`log_phi` for `family = "sltb"`, `log_sigma_e` for
#' `family = "gaussian"`). This is the numeric escape hatch consumed by tooling
#' such as `car::deltaMethod`.
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return Named numeric vector.
#' @export
coef.beezdiscounting_tmb <- function(object, ...) {
  object$model$coefficients
}

#' Extract fixed effects from a TMB discounting model
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return Named numeric vector (identical to `coef()`).
#' @export
fixef.beezdiscounting_tmb <- function(object, ...) {
  coef(object)
}

#' Extract subject-level random effects from a TMB discounting model
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return Data frame with `id`, the standardized random-intercept deviate
#'   `u_i` (such that `log k_i = X beta + sigma_u * u_i`), and the resolved
#'   per-subject discount rate `k`. A `phi` column is included only when the
#'   fit carries a subject-level precision (population-`phi` MVP omits it).
#' @export
ranef.beezdiscounting_tmb <- function(object, ...) {
  sp <- object$subject_pars
  keep <- intersect(c("id", "u_i", "k", "phi"), names(sp))
  out <- sp[, keep, drop = FALSE]
  rownames(out) <- NULL
  out
}
```

`fixef` and `ranef` are generics re-exported from `nlme`. Add to the package
roxygen (Phase F or DESCRIPTION should already `Imports: nlme`):

```r
#' @importFrom nlme fixef
#' @export
nlme::fixef

#' @importFrom nlme ranef
#' @export
nlme::ranef
```

Place these two re-export blocks in `R/beezdiscounting-package.R` (or any
collated R file) if not already present from Phase F.

- [ ] **Step 4: document + run (expect pass)**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: NAMESPACE gains `S3method(coef,beezdiscounting_tmb)`,
`S3method(fixef,beezdiscounting_tmb)`, `S3method(ranef,beezdiscounting_tmb)`,
and `export(fixef)` / `export(ranef)`; `coef / fixef` and `ranef` blocks pass.

- [ ] **Step 5: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R R/beezdiscounting-package.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): coef/fixef/ranef for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.3: `predict` (type = response/parameters; level = subject/population)

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("predict", {
  fit <- fit_for_methods()

  it("type = 'parameters' returns the subject_pars tibble", {
    pp <- predict(fit, type = "parameters")
    expect_s3_class(pp, "tbl_df")
    expect_true(all(c("id", "u_i", "k") %in% names(pp)))
    expect_equal(nrow(pp), fit$param_info$n_subjects)
  })

  it("type = 'response', level = 'subject' adds .fitted in [0,1]", {
    pr <- predict(fit, type = "response", level = "subject")
    expect_true(".fitted" %in% names(pr))
    expect_equal(nrow(pr), nrow(fit$data))
    expect_true(all(pr$.fitted >= 0 & pr$.fitted <= 1))
  })

  it("level = 'subject' equals the discounting fn at the subject's k", {
    pr <- predict(fit, type = "response", level = "subject")
    sp <- fit$subject_pars
    k_by_id <- stats::setNames(sp$k, sp$id)
    x <- fit$data[[fit$param_info$x_var]]
    id <- as.character(fit$data[[fit$param_info$id_var]])
    mu <- 1 / (1 + k_by_id[id] * x)            # mazur
    expect_equal(pr$.fitted, unname(mu), tolerance = 1e-8)
  })

  it("level = 'population' uses RE = 0 and needs no id column", {
    nd <- data.frame(x = c(7, 30, 180, 365))
    pr <- predict(fit, newdata = nd, type = "response", level = "population")
    expect_true("predict.fixed" %in% names(pr))
    k_pop <- exp(unname(fit$model$coefficients["beta_k"][1]))
    expect_equal(pr$predict.fixed, 1 / (1 + k_pop * nd$x), tolerance = 1e-8)
  })

  it("level = c('population','subject') returns both columns", {
    nd <- fit$data
    pr <- predict(fit, newdata = nd, type = "response",
                  level = c("population", "subject"))
    expect_true(all(c("predict.fixed", "predict.id") %in% names(pr)))
  })

  it("rejects a numeric nlme-style level", {
    expect_error(predict(fit, level = 1), "should be one of")
  })

  it("type = 'parameters' for an exponential fit yields .fitted = exp(-k x)", {
    set.seed(7)
    d2 <- .simulate_dd_ip_mixed(n_subjects = 25, family = "gaussian",
                                equation = "exponential", seed = 7)
    f2 <- fit_dd_tmb(d2, equation = "exponential", family = "gaussian",
                     random_effects = k ~ 1, verbose = 0)
    pr <- predict(f2, type = "response", level = "population",
                  newdata = data.frame(x = c(7, 365)))
    k_pop <- exp(unname(f2$model$coefficients["beta_k"][1]))
    expect_equal(pr$predict.fixed, exp(-k_pop * c(7, 365)), tolerance = 1e-8)
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — no `predict` method for class.

- [ ] **Step 3: implement (with the per-row mean back end)**

Append to `R/dd-tmb-methods.R`:

```r
# --- predict ---

#' Discounting function value mu = E[y] for a vector of k and x
#'
#' Mazur: `mu = 1 / (1 + k * x)`; exponential: `mu = exp(-k * x)`. Guards mu to
#' `[1e-6, 1 - 1e-6]` to match the C++ template.
#' @keywords internal
.dd_discount_mu <- function(k, x, equation) {
  mu <- switch(equation,
    mazur = 1 / (1 + k * x),
    exponential = exp(-k * x),
    stop("unknown equation '", equation, "'", call. = FALSE)
  )
  pmin(pmax(mu, 1e-6), 1 - 1e-6)
}

#' Rebuild the per-row log-k linear predictor from newdata
#'
#' Reconstructs the fixed-effect design from `newdata` using the stored RHS,
#' multiplies by `beta_k`, and (for `level = "subject"`) adds each subject's
#' `sigma_u * u_i` deviate looked up by id. Returns per-row `k`.
#'
#' @param object A `beezdiscounting_tmb` fit.
#' @param newdata Data frame with the x var, factor/covariate cols, and (for
#'   `level = "subject"`) the id var.
#' @param level `"subject"` or `"population"`.
#' @return Numeric vector of per-row `k`, length `nrow(newdata)`.
#' @keywords internal
.dd_tmb_predict_k <- function(object, newdata, level = c("subject", "population")) {
  level <- match.arg(level)
  pinfo <- object$param_info
  coefs <- object$model$coefficients
  beta_k <- unname(coefs[names(coefs) == "beta_k"])

  rhs <- object$formula_details$rhs %||% ~ 1
  Xnew <- stats::model.matrix(rhs, data = newdata)
  # Align columns to the fitted design (guards factor-level/order drift).
  fit_cn <- colnames(object$formula_details$X)
  miss <- setdiff(fit_cn, colnames(Xnew))
  if (length(miss) > 0L) {
    for (m in miss) Xnew <- cbind(Xnew, stats::setNames(rep(0, nrow(Xnew)), NULL))
    colnames(Xnew)[(ncol(Xnew) - length(miss) + 1L):ncol(Xnew)] <- miss
  }
  Xnew <- Xnew[, fit_cn, drop = FALSE]

  eta <- as.numeric(Xnew %*% beta_k)

  if (level == "subject") {
    id_var <- pinfo$id_var
    if (is.null(newdata[[id_var]])) {
      cli::cli_abort(c(
        "{.code level = \"subject\"} needs the id column {.val {id_var}} in {.arg newdata}.",
        i = "Use {.code level = \"population\"} for the random-effects-at-zero curve."
      ))
    }
    sp <- object$subject_pars
    sigma_u <- exp(unname(coefs[["log_sigma_u"]]))
    u_by_id <- stats::setNames(sp$u_i, as.character(sp$id))
    u_row <- u_by_id[as.character(newdata[[id_var]])]
    if (anyNA(u_row)) {
      cli::cli_abort("newdata contains ids not present in the fit.")
    }
    eta <- eta + sigma_u * unname(u_row)
  }
  exp(eta)
}

#' Predict from a TMB mixed-effects discounting model
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param newdata Optional data frame. `NULL` uses the fitting data.
#' @param type `"response"` (fitted indifference proportions on `[0,1]`) or
#'   `"parameters"` (the per-subject parameter tibble).
#' @param level For `type = "response"`: `"subject"` (default; conditions on
#'   each subject's random intercept, requires the id column) and/or
#'   `"population"` (random effects at zero — the population-mean curve, no id
#'   column needed). Pass `c("population", "subject")` for both. A numeric
#'   nlme-style level is rejected.
#' @param ... Unused.
#' @return For `type = "parameters"`, the subject-parameter tibble. For
#'   `type = "response"`, `newdata` plus a fitted column: `.fitted` when only
#'   `"subject"`; `predict.fixed` / `predict.id` when `"population"` is
#'   requested (nlme-style names so nlme-based plotting runs unchanged).
#' @export
predict.beezdiscounting_tmb <- function(object,
                                        newdata = NULL,
                                        type = c("response", "parameters"),
                                        level = "subject",
                                        ...) {
  type <- match.arg(type)
  level <- match.arg(level, c("subject", "population"), several.ok = TRUE)

  if (type == "parameters") {
    return(tibble::as_tibble(object$subject_pars))
  }

  equation <- object$param_info$equation
  x_var <- object$param_info$x_var
  if (is.null(newdata)) newdata <- object$data
  out <- tibble::as_tibble(newdata)
  x <- newdata[[x_var]]

  if (identical(level, "subject")) {
    k_row <- .dd_tmb_predict_k(object, newdata, level = "subject")
    out$.fitted <- .dd_discount_mu(k_row, x, equation)
    return(out)
  }
  if ("population" %in% level) {
    k_pop <- .dd_tmb_predict_k(object, newdata, level = "population")
    out$predict.fixed <- .dd_discount_mu(k_pop, x, equation)
  }
  if ("subject" %in% level) {
    k_sub <- .dd_tmb_predict_k(object, newdata, level = "subject")
    out$predict.id <- .dd_discount_mu(k_sub, x, equation)
  }
  out
}
```

- [ ] **Step 4: document + run (expect pass)**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: `predict` block passes (8 `it`s); NAMESPACE gains
`S3method(predict,beezdiscounting_tmb)`.

- [ ] **Step 5: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): predict (response/parameters, subject/population) for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.4: `augment` (.fitted / .resid / .std_resid; supports newdata)

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("augment", {
  fit <- fit_for_methods()

  it("adds .fitted, .resid, .std_resid and preserves rows", {
    aug <- augment(fit)
    expect_s3_class(aug, "tbl_df")
    expect_true(all(c(".fitted", ".resid", ".std_resid") %in% names(aug)))
    expect_equal(nrow(aug), nrow(fit$data))
    y <- fit$data[[fit$param_info$y_var]]
    expect_equal(aug$.resid, y - aug$.fitted, tolerance = 1e-10)
  })

  it("std_resid scales resid by the SLT delta-method SD per row", {
    aug <- augment(fit)
    # .std_resid = .resid / sd_i, sd_i finite and positive
    expect_true(all(is.finite(aug$.std_resid)))
    ratio <- aug$.resid / aug$.std_resid
    expect_true(all(ratio > 0))
  })

  it("supports newdata and computes resid against the y_var in newdata", {
    nd <- fit$data[1:10, , drop = FALSE]
    aug <- augment(fit, newdata = nd)
    expect_equal(nrow(aug), 10L)
    expect_true(all(c(".fitted", ".resid") %in% names(aug)))
  })

  it("gaussian fit uses constant sigma_e for .std_resid", {
    set.seed(9)
    d2 <- .simulate_dd_ip_mixed(n_subjects = 25, family = "gaussian",
                                equation = "mazur", seed = 9)
    f2 <- fit_dd_tmb(d2, equation = "mazur", family = "gaussian",
                     random_effects = k ~ 1, verbose = 0)
    aug <- augment(f2)
    sigma_e <- exp(f2$model$coefficients[["log_sigma_e"]])
    expect_equal(aug$.std_resid, aug$.resid / sigma_e, tolerance = 1e-10)
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — no `augment` method for class.

- [ ] **Step 3: implement the shared fitted/resid back end + augment**

Append to `R/dd-tmb-methods.R`:

```r
# --- fitted / residuals / augment ---

#' Per-row response SD on the [0,1] scale for standardized residuals
#'
#' Gaussian: constant `sigma_e = exp(log_sigma_e)`. SLT-beta: the
#' delta-method SLT SD `s * sqrt(mu (1 - mu) / (phi + 1))` (the verified SLT
#' moment from the spec, keeping the `s^2` factor), so residuals near the
#' bounds are down-weighted as intended.
#' @keywords internal
.dd_tmb_response_sd <- function(object, mu) {
  coefs <- object$model$coefficients
  family <- object$param_info$family
  if (family == "gaussian") {
    return(rep(exp(coefs[["log_sigma_e"]]), length(mu)))
  }
  s <- 1.0000001
  phi <- exp(coefs[["log_phi"]])
  s * sqrt(mu * (1 - mu) / (phi + 1))
}

#' Fitted values + residuals shared by fitted(), residuals(), augment()
#' @keywords internal
.dd_tmb_fitted_resid <- function(object, newdata = NULL,
                                 level = c("subject", "population")) {
  level <- match.arg(level)
  data_used <- if (is.null(newdata)) object$data else newdata
  pred <- predict(object, newdata = data_used, type = "response", level = level)
  fitted_vals <- if (level == "population") pred$predict.fixed else pred$.fitted
  y_obs <- data_used[[object$param_info$y_var]]
  list(.fitted = fitted_vals, .resid = y_obs - fitted_vals, data = data_used)
}

#' Fitted values for a beezdiscounting_tmb fit
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param level `"subject"` (default) or `"population"`.
#' @param ... Unused.
#' @return Numeric vector of length `nobs(object)`.
#' @export
fitted.beezdiscounting_tmb <- function(object,
                                       level = c("subject", "population"),
                                       ...) {
  level <- match.arg(level)
  .dd_tmb_fitted_resid(object, level = level)$.fitted
}

#' Residuals for a beezdiscounting_tmb fit
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param type `"response"` (default) or `"pearson"` (divided by the per-row
#'   response SD; see [augment.beezdiscounting_tmb()]).
#' @param level `"subject"` (default) or `"population"`.
#' @param ... Unused.
#' @return Numeric vector of length `nobs(object)`.
#' @export
residuals.beezdiscounting_tmb <- function(object,
                                          type = c("response", "pearson"),
                                          level = c("subject", "population"),
                                          ...) {
  type <- match.arg(type)
  level <- match.arg(level)
  fr <- .dd_tmb_fitted_resid(object, level = level)
  if (type == "response") return(fr$.resid)
  fr$.resid / .dd_tmb_response_sd(object, fr$.fitted)
}

#' Augment a beezdiscounting_tmb model
#'
#' @param x A `beezdiscounting_tmb` object.
#' @param newdata Optional data frame. `NULL` uses the fitting data.
#' @param ... Unused.
#' @return A tibble of the data plus `.fitted` (subject-conditional fitted
#'   indifference proportion), `.resid` (`y - .fitted` on the `[0,1]`
#'   response scale), and `.std_resid` (`.resid` divided by the per-row
#'   response SD: constant `sigma_e` for `gaussian`, the SLT delta-method SD
#'   for `sltb`).
#' @export
augment.beezdiscounting_tmb <- function(x, newdata = NULL, ...) {
  fr <- .dd_tmb_fitted_resid(x, newdata = newdata, level = "subject")
  out <- tibble::as_tibble(fr$data)
  out$.fitted <- fr$.fitted
  out$.resid <- fr$.resid
  out$.std_resid <- fr$.resid / .dd_tmb_response_sd(x, fr$.fitted)
  out
}
```

`augment` / `tidy` / `glance` are generics from `generics` (re-exported by
broom). Add re-exports (place in `R/beezdiscounting-package.R` if Phase P/F
has not already added them):

```r
#' @importFrom generics augment
#' @export
generics::augment

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance
```

- [ ] **Step 4: document + run (expect pass)**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: `augment` block passes (4 `it`s); NAMESPACE gains
`S3method(augment,beezdiscounting_tmb)`,
`S3method(fitted,beezdiscounting_tmb)`,
`S3method(residuals,beezdiscounting_tmb)`, and `export(augment)` /
`export(tidy)` / `export(glance)`.

- [ ] **Step 5: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R R/beezdiscounting-package.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): augment/fitted/residuals for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.5: `tidy` (report_space natural/log10/internal; broom contract)

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

Depends on Phase P's `.dd_transform_coef_table()` (the beezdemand port with
regex retargeted to `^k($|_)|^s($|_)|^phi($|_)`) and
`.dd_transform_est_se()`.

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("tidy", {
  fit <- fit_for_methods()

  it("returns the broom column contract", {
    td <- tidy(fit)
    expect_s3_class(td, "tbl_df")
    expect_true(all(c("term", "estimate", "std.error", "statistic",
                      "p.value", "component", "estimate_scale",
                      "term_display") %in% names(td)))
  })

  it("fixed-effect rows carry component = 'fixed', variance rows 'variance'", {
    td <- tidy(fit)
    expect_true(any(td$component == "fixed"))
    expect_true(any(td$component == "variance"))
    k_row <- td[td$term == "k:(Intercept)", ]
    expect_equal(nrow(k_row), 1L)
    expect_equal(k_row$component, "fixed")
  })

  it("internal space leaves the log-k estimate untransformed", {
    td <- tidy(fit, report_space = "internal")
    k_row <- td[td$term == "k:(Intercept)", ]
    expect_equal(k_row$estimate, unname(fit$model$coefficients["beta_k"][1]),
                 tolerance = 1e-10)
    expect_equal(k_row$estimate_scale, "log")
  })

  it("natural space exponentiates the intercept to k and rescales SE", {
    td_int <- tidy(fit, report_space = "internal")
    td_nat <- tidy(fit, report_space = "natural")
    ki <- td_int[td_int$term == "k:(Intercept)", ]
    kn <- td_nat[td_nat$term == "k:(Intercept)", ]
    expect_equal(kn$estimate, exp(ki$estimate), tolerance = 1e-8)
    expect_equal(kn$std.error, exp(ki$estimate) * ki$std.error, tolerance = 1e-8)
    expect_equal(kn$estimate_scale, "natural")
  })

  it("log10 space divides the log estimate by log(10)", {
    td_int <- tidy(fit, report_space = "internal")
    td_l10 <- tidy(fit, report_space = "log10")
    ki <- td_int[td_int$term == "k:(Intercept)", ]
    kl <- td_l10[td_l10$term == "k:(Intercept)", ]
    expect_equal(kl$estimate, ki$estimate / log(10), tolerance = 1e-8)
    expect_equal(kl$estimate_scale, "log10")
  })

  it("statistic/p.value are estimation-scale invariant to report_space", {
    a <- tidy(fit, report_space = "internal")
    b <- tidy(fit, report_space = "natural")
    ka <- a[a$term == "k:(Intercept)", ]
    kb <- b[b$term == "k:(Intercept)", ]
    expect_equal(ka$statistic, kb$statistic, tolerance = 1e-10)
    expect_equal(ka$p.value, kb$p.value, tolerance = 1e-10)
  })

  it("effects = 'fixed' drops the variance rows", {
    td <- tidy(fit, effects = "fixed")
    expect_true(all(td$component == "fixed"))
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — no `tidy` method for class.

- [ ] **Step 3: implement**

Append to `R/dd-tmb-methods.R`:

```r
# --- tidy / glance ---

#' Tidy a beezdiscounting_tmb model
#'
#' @param x A `beezdiscounting_tmb` object.
#' @param effects `"fixed"` (log-k fixed-effect rows), `"ran_pars"` (the RE SD
#'   and the auxiliary precision/scale), or both (default).
#' @param report_space `"natural"`, `"log10"`, or `"internal"` — scale for the
#'   fixed-effect rows. `estimate`/`std.error` follow this scale;
#'   `statistic`/`p.value` are always on the estimation (log-k) scale.
#'   Variance-component rows are unaffected.
#' @param ... Unused.
#' @return A tibble with `term`, `estimate`, `std.error`, `statistic`,
#'   `p.value`, `component`, `estimate_scale`, `term_display`. Fixed-effect
#'   rows carry `component == "fixed"`; variance rows `component == "variance"`.
#' @export
tidy.beezdiscounting_tmb <- function(x,
                                     effects = c("fixed", "ran_pars"),
                                     report_space = c("natural", "log10", "internal"),
                                     ...) {
  effects <- match.arg(effects, several.ok = TRUE)
  report_space <- match.arg(report_space)

  result <- tibble::tibble()

  if ("fixed" %in% effects) {
    coefs <- x$model$coefficients
    se <- .dd_tmb_model_se(x)
    nms <- names(coefs)
    tn <- .dd_tmb_build_term_names(x, nms)

    component <- ifelse(nms == "beta_k", "fixed", "variance")
    estimate_scale <- rep("log", length(nms))
    # sigma_u/phi/sigma_e are reported on the log scale internally too; the
    # fixed rows we keep are beta_k only, which are log-k.

    z_val <- coefs / se
    p_val <- 2 * stats::pnorm(-abs(z_val))

    fixed <- tibble::tibble(
      term = tn$term,
      estimate = unname(coefs),
      std.error = unname(se),
      statistic = unname(z_val),
      p.value = unname(p_val),
      component = component,
      estimate_scale = estimate_scale,
      term_display = tn$term
    )
    fixed <- fixed[fixed$component == "fixed", , drop = FALSE]

    fixed <- .dd_transform_coef_table(
      coef_tbl = fixed,
      report_space = report_space,
      internal_space = "log"
    )

    result <- dplyr::bind_rows(result, fixed)
  }

  if ("ran_pars" %in% effects) {
    vc <- .dd_tmb_variance_components(x)
    ran <- tibble::tibble(
      term = vc$Component,
      estimate = vc$Estimate,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      component = "variance",
      estimate_scale = vc$Scale,
      term_display = vc$Component
    )
    result <- dplyr::bind_rows(result, ran)
  }

  if (isFALSE(x$converged) || isFALSE(x$se_available)) {
    attr(result, "se_warning") <-
      "Standard errors / convergence are unreliable; CIs and p-values may be invalid."
  }

  result
}
```

Also add the variance-component formatter (used by `tidy(effects="ran_pars")`,
`summary()`, and `confint()`):

```r
#' Random-effect + auxiliary variance components for a TMB discounting fit
#'
#' Reports the random-intercept SD on log-k (`sigma_u = exp(log_sigma_u)`,
#' on the natural-log scale, divided by log(10) for log10 comparability with
#' nlme::VarCorr), and the auxiliary parameter: `phi` (sltb) or `sigma_e`
#' (gaussian), both back-transformed to their natural scale.
#' @keywords internal
.dd_tmb_variance_components <- function(object) {
  coefs <- object$model$coefficients
  family <- object$param_info$family
  ln10 <- log(10)

  rows <- list(
    data.frame(
      Component = "sigma_u (log10-k RE SD)",
      Estimate = exp(coefs[["log_sigma_u"]]) / ln10,
      Scale = "log10",
      stringsAsFactors = FALSE
    )
  )
  if (family == "sltb") {
    rows[[length(rows) + 1L]] <- data.frame(
      Component = "phi (precision)",
      Estimate = exp(coefs[["log_phi"]]),
      Scale = "natural",
      stringsAsFactors = FALSE
    )
  } else {
    rows[[length(rows) + 1L]] <- data.frame(
      Component = "sigma_e (Residual SD)",
      Estimate = exp(coefs[["log_sigma_e"]]),
      Scale = "natural",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
```

- [ ] **Step 4: document + run (expect pass)**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: `tidy` block passes (7 `it`s); NAMESPACE gains
`S3method(tidy,beezdiscounting_tmb)`.

- [ ] **Step 5: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): tidy + variance-components for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.6: `glance` (1-row model summary)

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("glance", {
  fit <- fit_for_methods()

  it("returns one row with the canonical columns and backend string", {
    g <- glance(fit)
    expect_s3_class(g, "tbl_df")
    expect_equal(nrow(g), 1L)
    expect_named(g, c("model_class", "backend", "equation", "family",
                      "nobs", "n_subjects", "n_random_effects",
                      "converged", "logLik", "AIC", "BIC"))
    expect_equal(g$model_class, "beezdiscounting_tmb")
    expect_equal(g$backend, "TMB_mixed")
    expect_equal(g$equation, "mazur")
    expect_equal(g$family, "sltb")
    expect_equal(g$nobs, fit$param_info$n_obs)
    expect_equal(g$n_subjects, fit$param_info$n_subjects)
    expect_equal(g$logLik, fit$loglik)
    expect_equal(g$AIC, fit$AIC)
    expect_equal(g$BIC, fit$BIC)
  })

  it("reports family = gaussian for a gaussian fit", {
    set.seed(11)
    d2 <- .simulate_dd_ip_mixed(n_subjects = 25, family = "gaussian",
                                equation = "mazur", seed = 11)
    f2 <- fit_dd_tmb(d2, equation = "mazur", family = "gaussian",
                     random_effects = k ~ 1, verbose = 0)
    expect_equal(glance(f2)$family, "gaussian")
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — no `glance` method for class.

- [ ] **Step 3: implement**

Append to `R/dd-tmb-methods.R`:

```r
#' Glance at a beezdiscounting_tmb model
#'
#' @param x A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return A one-row tibble: `model_class` (`"beezdiscounting_tmb"`),
#'   `backend` (`"TMB_mixed"`), `equation`, `family`, `nobs`, `n_subjects`,
#'   `n_random_effects`, `converged`, `logLik`, `AIC`, `BIC`.
#' @export
glance.beezdiscounting_tmb <- function(x, ...) {
  tibble::tibble(
    model_class = "beezdiscounting_tmb",
    backend = "TMB_mixed",
    equation = x$param_info$equation,
    family = x$param_info$family,
    nobs = x$param_info$n_obs,
    n_subjects = x$param_info$n_subjects,
    n_random_effects = x$param_info$n_random_effects,
    converged = x$converged,
    logLik = x$loglik,
    AIC = x$AIC,
    BIC = x$BIC
  )
}
```

- [ ] **Step 4: document + run (expect pass)**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: `glance` block passes; NAMESPACE gains
`S3method(glance,beezdiscounting_tmb)`.

- [ ] **Step 5: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): glance for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.7: `confint` (Wald; report_space natural/internal; parm filter)

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("confint", {
  fit <- fit_for_methods()

  it("returns term/estimate/conf.low/conf.high/level with internal default", {
    ci <- confint(fit)
    expect_s3_class(ci, "tbl_df")
    expect_named(ci, c("term", "estimate", "conf.low", "conf.high", "level"))
    expect_true(all(ci$conf.low <= ci$estimate))
    expect_true(all(ci$estimate <= ci$conf.high))
    expect_true(all(ci$level == 0.95))
  })

  it("Wald interval = estimate +/- z*se on the internal scale", {
    ci <- confint(fit, parm = "k:(Intercept)")
    co <- fit$model$coefficients["beta_k"][1]
    se <- .dd_tmb_model_se(fit)["beta_k"]
    z <- stats::qnorm(0.975)
    expect_equal(ci$estimate, unname(co), tolerance = 1e-10)
    expect_equal(ci$conf.low, unname(co - z * se), tolerance = 1e-8)
    expect_equal(ci$conf.high, unname(co + z * se), tolerance = 1e-8)
  })

  it("report_space = 'natural' exponentiates beta_k bounds", {
    ci_i <- confint(fit, parm = "k:(Intercept)", report_space = "internal")
    ci_n <- confint(fit, parm = "k:(Intercept)", report_space = "natural")
    expect_equal(ci_n$estimate, exp(ci_i$estimate), tolerance = 1e-8)
    expect_equal(ci_n$conf.low, exp(ci_i$conf.low), tolerance = 1e-8)
    expect_equal(ci_n$conf.high, exp(ci_i$conf.high), tolerance = 1e-8)
  })

  it("parm filters by display name OR raw name", {
    ci_disp <- confint(fit, parm = "k:(Intercept)")
    ci_raw <- confint(fit, parm = "beta_k")
    expect_equal(nrow(ci_disp), 1L)
    expect_equal(nrow(ci_raw), 1L)
    expect_equal(ci_disp$estimate, ci_raw$estimate, tolerance = 1e-12)
  })

  it("honours a non-default confidence level", {
    ci90 <- confint(fit, parm = "k:(Intercept)", level = 0.90)
    ci95 <- confint(fit, parm = "k:(Intercept)", level = 0.95)
    expect_true(ci90$conf.low > ci95$conf.low)
    expect_true(ci90$conf.high < ci95$conf.high)
    expect_true(all(ci90$level == 0.90))
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — no `confint` method for class.

- [ ] **Step 3: implement**

Append to `R/dd-tmb-methods.R`:

```r
# --- confint ---

#' Confidence intervals for a TMB discounting model
#'
#' Wald (Hessian-based) intervals `estimate +/- z * se`.
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param parm Optional character vector filtering by display name
#'   (`"k:(Intercept)"`) or raw name (`"beta_k"`, `"log_sigma_u"`,
#'   `"log_phi"`, `"log_sigma_e"`). `NULL` returns all coefficients.
#' @param level Confidence level (default `0.95`).
#' @param report_space `"internal"` (default; all coefficients on their
#'   estimation/log scale) or `"natural"` (exponentiate `beta_k` so the
#'   intercept is `k` at the reference level and non-intercept terms are
#'   multiplicative fold-changes; variance/aux params stay internal).
#' @param ... Unused.
#' @return A tibble with `term`, `estimate`, `conf.low`, `conf.high`, `level`.
#' @export
confint.beezdiscounting_tmb <- function(object,
                                        parm = NULL,
                                        level = 0.95,
                                        report_space = c("internal", "natural"),
                                        ...) {
  report_space <- match.arg(report_space)
  coefs <- object$model$coefficients
  se_vec <- .dd_tmb_model_se(object)
  nms <- names(coefs)
  tn <- .dd_tmb_build_term_names(object, nms)
  term <- tn$term

  if (!is.null(parm)) {
    keep <- term %in% parm | nms %in% parm
    coefs <- coefs[keep]; se_vec <- se_vec[keep]
    nms <- nms[keep]; term <- term[keep]
  }

  estimates <- coefs
  z <- stats::qnorm((1 + level) / 2)
  conf_low <- coefs - z * se_vec
  conf_high <- coefs + z * se_vec

  if (report_space == "natural") {
    k_idx <- which(nms == "beta_k")
    if (length(k_idx) > 0L) {
      estimates[k_idx] <- exp(coefs[k_idx])
      conf_low[k_idx] <- exp(conf_low[k_idx])
      conf_high[k_idx] <- exp(conf_high[k_idx])
    }
  }

  tibble::tibble(
    term = term,
    estimate = unname(estimates),
    conf.low = unname(conf_low),
    conf.high = unname(conf_high),
    level = level
  )
}
```

- [ ] **Step 4: document + run (expect pass)**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: `confint` block passes (5 `it`s); NAMESPACE gains
`S3method(confint,beezdiscounting_tmb)`.

- [ ] **Step 5: commit**

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): confint (Wald, natural/internal) for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task M.8: `summary` / `print` / `print.summary`

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-tmb-methods.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R`

- [ ] **Step 1: write failing tests**

Append to `tests/testthat/test-dd-tmb-methods.R`:

```r
describe("summary / print", {
  fit <- fit_for_methods()

  it("summary() returns the class with coefficients + variance components", {
    s <- summary(fit)
    expect_s3_class(s, "summary.beezdiscounting_tmb")
    expect_equal(s$backend, "TMB_mixed")
    expect_equal(s$equation, "mazur")
    expect_equal(s$family, "sltb")
    expect_s3_class(s$coefficients, "tbl_df")
    expect_true(all(c("term", "estimate", "std.error", "statistic",
                      "p.value") %in% names(s$coefficients)))
    expect_true(!is.null(s$variance_components))
    expect_equal(s$n_subjects, fit$param_info$n_subjects)
  })

  it("summary() honours report_space for the fixed-effect estimates", {
    s_nat <- summary(fit, report_space = "natural")
    krow <- s_nat$coefficients[s_nat$coefficients$term == "k:(Intercept)", ]
    expect_equal(krow$estimate,
                 exp(unname(fit$model$coefficients["beta_k"][1])),
                 tolerance = 1e-8)
  })

  it("print() and print.summary() run without error and return invisibly", {
    expect_invisible(print(fit))
    expect_output(print(fit), "TMB Mixed-Effects Discounting Model")
    s <- summary(fit)
    expect_invisible(print(s))
    expect_output(print(s), "TMB_mixed")
  })

  it("summary notes flag non-convergence / missing SEs", {
    s <- summary(fit)
    expect_type(s$notes, "character")
  })
})
```

- [ ] **Step 2: run (expect fail)**

```
Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: failures — no `summary` method for class.

- [ ] **Step 3: implement summary + both print methods**

Append to `R/dd-tmb-methods.R`:

```r
# --- summary / print ---

#' Summarize a TMB mixed-effects discounting fit
#'
#' @param object A `beezdiscounting_tmb` object.
#' @param report_space Scale for the fixed-effect (`beta_k`) estimates /
#'   std.errors: `"natural"` (default), `"log10"`, or `"internal"`.
#'   `statistic`/`p.value` are always on the estimation (log-k) scale.
#' @param ... Unused.
#' @return An object of class `summary.beezdiscounting_tmb`.
#' @export
summary.beezdiscounting_tmb <- function(object,
                                        report_space = c("natural", "log10", "internal"),
                                        ...) {
  report_space <- match.arg(report_space)

  coefs <- object$model$coefficients
  se_vec <- .dd_tmb_model_se(object)
  nms <- names(coefs)
  tn <- .dd_tmb_build_term_names(object, nms)

  component <- ifelse(nms == "beta_k", "fixed", "variance")
  estimate_scale <- rep("log", length(nms))

  z_val <- coefs / se_vec
  p_val <- 2 * stats::pnorm(-abs(z_val))

  coefficients <- tibble::tibble(
    term = tn$term,
    estimate = unname(coefs),
    std.error = unname(se_vec),
    statistic = unname(z_val),
    p.value = unname(p_val),
    component = component,
    estimate_scale = estimate_scale,
    term_display = tn$term
  )
  coefficients <- coefficients[coefficients$component == "fixed", , drop = FALSE]
  coefficients <- .dd_transform_coef_table(
    coef_tbl = coefficients,
    report_space = report_space,
    internal_space = "log"
  )

  vc <- .dd_tmb_variance_components(object)

  notes <- character(0)
  if (!object$converged) {
    notes <- c(notes, "WARNING: Model did not converge.")
  }
  if (isFALSE(object$se_available)) {
    notes <- c(notes, "Standard errors unavailable (sdreport failed); CIs will be NA.")
  }
  if (length(object$opt_warnings %||% character(0)) > 0) {
    notes <- c(notes, sprintf("Optimizer produced %d warning(s) during fitting.",
                              length(object$opt_warnings)))
  }
  if (!is.null(object$param_info$factors) && length(object$param_info$factors) > 0) {
    notes <- c(notes,
      "Population k reflects the reference level. Use get_dd_param_emms() for per-group estimates.")
  }

  structure(
    list(
      call = object$opt$message,
      model_class = "beezdiscounting_tmb",
      backend = "TMB_mixed",
      equation = object$param_info$equation,
      family = object$param_info$family,
      coefficients = coefficients,
      variance_components = vc,
      n_subjects = object$param_info$n_subjects,
      nobs = object$param_info$n_obs,
      converged = object$converged,
      logLik = object$loglik,
      AIC = object$AIC,
      BIC = object$BIC,
      notes = notes
    ),
    class = "summary.beezdiscounting_tmb"
  )
}

#' Print a TMB mixed-effects discounting fit
#'
#' @param x A `beezdiscounting_tmb` object.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.beezdiscounting_tmb <- function(x, ...) {
  cat("\nTMB Mixed-Effects Discounting Model\n\n")
  cat("Equation:", x$param_info$equation, "\n")
  cat("Family:", x$param_info$family, "\n")
  cat("Convergence:", ifelse(x$converged, "Yes", "No"), "\n")
  cat("Number of subjects:", x$param_info$n_subjects, "\n")
  cat("Number of observations:", x$param_info$n_obs, "\n")
  cat("Random effects:", x$param_info$n_random_effects, "(k ~ 1)\n")
  cat("Log-likelihood:", round(x$loglik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")
  cat("\nFixed Effects (log k):\n")
  co <- x$model$coefficients
  tn <- .dd_tmb_build_term_names(x)
  fe <- co[tn$k_idx]
  names(fe) <- tn$term[tn$k_idx]
  print(round(fe, 4))
  cat("\nUse summary() for full results.\n")
  invisible(x)
}

#' Print a TMB discounting model summary
#'
#' @param x A `summary.beezdiscounting_tmb` object.
#' @param digits Significant digits.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.summary.beezdiscounting_tmb <- function(x, digits = 4, ...) {
  cat("\nTMB Mixed-Effects Discounting Model Summary\n")
  cat(strrep("=", 50), "\n\n")
  cat("Equation:", x$equation, "\n")
  cat("Family:", x$family, "\n")
  cat("Backend:", x$backend, "\n")
  cat("Convergence:", ifelse(x$converged, "Yes", "No"), "\n")
  cat("Subjects:", x$n_subjects, " Observations:", x$nobs, "\n\n")

  cat("--- Fixed Effects (log k) ---\n")
  cd <- as.data.frame(x$coefficients[, c("term", "estimate", "std.error",
                                         "statistic", "p.value")])
  cd$estimate <- round(cd$estimate, digits)
  cd$std.error <- round(cd$std.error, digits)
  cd$statistic <- round(cd$statistic, digits)
  cd$p.value <- format.pval(cd$p.value, digits = 3)
  print(cd, row.names = FALSE)

  cat("\n--- Variance Components ---\n")
  vc <- x$variance_components
  vc$Estimate <- round(vc$Estimate, digits)
  print(vc, row.names = FALSE)

  cat("\n--- Fit Statistics ---\n")
  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "  BIC:", round(x$BIC, 2), "\n")

  if (length(x$notes) > 0) {
    cat("\nNotes:\n")
    for (note in x$notes) cat("  *", note, "\n")
  }
  invisible(x)
}
```

- [ ] **Step 4: document + run (expect pass)**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: `summary / print` block passes; NAMESPACE gains
`S3method(summary,beezdiscounting_tmb)`,
`S3method(print,beezdiscounting_tmb)`,
`S3method(print,summary.beezdiscounting_tmb)`.

- [ ] **Step 5: full-file run + commit**

```
Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting"); devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-tmb-methods.R")'
```
Expected: every `describe` block in the file passes; 0 failures, 0 warnings.

```
git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-tmb-methods.R tests/testthat/test-dd-tmb-methods.R NAMESPACE
git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(methods): summary/print/print.summary for beezdiscounting_tmb

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
## Phase E — emmeans & comparisons

Phase E builds the "EMMs without emmeans" averaging-matrix surface for the
single-parameter discounting model. Everything is computed for `k` on the
`log k` scale (linear in `beta_k`), then back-transformed by `exp()`. Because
the family (sltb/gaussian) only changes the C++ density and the auxiliary
parameter, this entire phase is **family-agnostic**: it touches only
`fit$model$coefficients[beta_k]`, `fit$sdr$cov.fixed`, `fit$formula_details$X`,
and `fit$param_info` (factors / continuous_covariates / factor_interaction /
data).

This is a structural simplification of the beezdemand TMB emmeans surface
(`/Users/brent/Dropbox/GIT/beezdemand/R/tmb-methods.R`): beezdemand carries two
parameters (`Q0`, `alpha`) with two design matrices (`X_q0`, `X_alpha`) and a
`param` switch; the discounting model has exactly one parameter (`k`) with one
design matrix (`X`), so the `param` argument and the per-parameter loop
collapse away. The averaging-matrix math (full-factorial grid →
`contrasts.arg = attr(fitted_X, "contrasts")` → reorder to fitted cols →
averaging matrix `A` → `ref_X = A %*% X_full`) is **identical** and ported
verbatim.

Preconditions from earlier phases (PINNED CONTRACT):
- `fit$formula_details$X` is the fitted log-k FE design matrix (from
  `.dd_tmb_build_design()`), with `attr(X, "contrasts")` set.
- `fit$param_info` has `factors`, `factor_interaction`, `continuous_covariates`,
  `id_var`, `x_var`, `y_var`.
- `fit$model$coefficients` includes a `beta_k` block (one entry per `ncol(X)`).
- `fit$opt$par` is the named optimizer vector (`beta_k` repeated ncol(X) times,
  then `log_sigma_u`, `log_aux`).
- `fit$sdr` is the `TMB::sdreport`; `fit$sdr$cov.fixed` is the fixed-effect
  covariance with the same `names(fit$opt$par)` ordering.
- `R/utils.R::build_fixed_rhs(factors, factor_interaction, continuous_covariates, data)`
  exists (ported from beezdemand `R/utils.R:657`) and returns the RHS string.
- `tibble`, `dplyr`, `cli`, `generics` (for the `tidy` generic) are in
  `Imports`.

---

### Task E.1: `.dd_build_emm_ref_grid()` — conditioned reference grid + averaging matrix

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-comparisons.R`
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R`

- [ ] **Step 1: Write the failing test for the ref-grid builder (non-circular: A and ref_X recomputed by hand).**

```r
# tests/testthat/test-dd-comparisons.R
skip_on_cran()
skip_if_not_installed("TMB")

describe(".dd_build_emm_ref_grid()", {

  it("returns is_intercept_only for a factor-free, covariate-free fit", {
    set.seed(1)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 40, family = "gaussian", equation = "mazur", seed = 1
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      multi_start = FALSE, verbose = 0)
    g <- .dd_build_emm_ref_grid(fit, at = NULL, factors_in_emm = NULL)
    expect_true(g$is_intercept_only)
    expect_null(g$ref_X)
    expect_null(g$level_combos)
  })

  it("builds ref_X = A %*% X_full with equal averaging weights over fitted cols", {
    set.seed(2)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.8, seed = 2
    )
    # n_conditions = 2 attaches a `condition` factor (PINNED CONTRACT simulator).
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)

    g <- .dd_build_emm_ref_grid(fit, at = NULL, factors_in_emm = NULL)
    expect_false(g$is_intercept_only)
    # Two observed cells (the two conditions), each a single full-grid row.
    expect_equal(nrow(g$ref_X), 2L)
    expect_equal(ncol(g$ref_X), ncol(fit$formula_details$X))
    expect_setequal(colnames(g$ref_X), colnames(fit$formula_details$X))

    # Non-circular check: rebuild X_full ourselves and confirm ref_X selects
    # the per-cell design rows (A = identity over observed cells when nothing
    # is marginalized).
    full <- expand.grid(condition = levels(fit$data$condition),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    full$condition <- factor(full$condition, levels = levels(fit$data$condition))
    X_full <- stats::model.matrix(~ condition, data = full,
                                  contrasts.arg = attr(fit$formula_details$X, "contrasts"))
    X_full <- X_full[, colnames(fit$formula_details$X), drop = FALSE]
    expect_equal(unname(as.matrix(g$ref_X)), unname(X_full), tolerance = 1e-12)
  })

  it("marginalizes an omitted factor with equal weights (A row sums to 1)", {
    set.seed(3)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.8, seed = 3
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    # factors_in_emm = character(0) marginalizes everything to one grand-mean row.
    g <- .dd_build_emm_ref_grid(fit, at = NULL, factors_in_emm = character(0))
    expect_equal(nrow(g$ref_X), 1L)
    # The grand-mean row is the equal-weight average of the per-condition rows.
    full <- expand.grid(condition = levels(fit$data$condition),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    full$condition <- factor(full$condition, levels = levels(fit$data$condition))
    X_full <- stats::model.matrix(~ condition, data = full,
                                  contrasts.arg = attr(fit$formula_details$X, "contrasts"))
    X_full <- X_full[, colnames(fit$formula_details$X), drop = FALSE]
    expect_equal(unname(as.numeric(g$ref_X[1, ])),
                 unname(colMeans(X_full)), tolerance = 1e-12)
  })
})
```

- [ ] **Step 2: Run the test, expect failure (function undefined).**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: errors `could not find function ".dd_build_emm_ref_grid"` (and `fit_dd_tmb`/`.simulate_dd_ip_mixed` already exist from earlier phases).

- [ ] **Step 3: Implement `.dd_build_emm_ref_grid()` and `.dd_resolve_retained_factors()` in `R/dd-comparisons.R`.**

  This is the single-design-matrix simplification of beezdemand
  `R/tmb-methods.R:3164-3398` (`.tmb_resolve_retained_factors` +
  `.tmb_build_emm_ref_grid`). Adaptations: drop the `param` argument and the
  `param == "Q0"`/`"alpha"` switch; the fitted design is always
  `fit_obj$formula_details$X`; fitted factors are always
  `fit_obj$param_info$factors`. The `.tmb_validate_at` machinery (collapse
  aliasing, per-param scope) is **omitted** — discounting has one parameter, so
  a plain at-validation suffices.

```r
# R/dd-comparisons.R

# Resolve user-requested retained factors against the fit's fitted factor set.
# Single-parameter analogue of beezdemand .tmb_resolve_retained_factors()
# (no collapse_levels / per-param suffixing in the discounting model).
.dd_resolve_retained_factors <- function(requested, fitted_factors) {
  unresolved <- setdiff(requested, fitted_factors)
  if (length(unresolved) > 0L) {
    cli::cli_abort(c(
      "{cli::qty(unresolved)}Requested factor{?s} {.val {unresolved}} {?is/are} not in the model.",
      "i" = "Fitted factors: {.val {fitted_factors}}."
    ))
  }
  unique(intersect(requested, fitted_factors))
}

# Validate the `at` list for the discounting EMM/comparison helpers. Aborts on
# unnamed entries, names not in (factors u continuous_covariates), factor values
# not observed, or non-finite continuous values; warns once on multi-value
# continuous entries (first value used). Single-parameter simplification of
# beezdemand .tmb_validate_at() (no param_scope / collapse alias logic).
.dd_validate_at <- function(fit_obj, at) {
  if (is.null(at)) return(invisible(NULL))
  if (is.null(names(at)) || any(!nzchar(names(at)))) {
    cli::cli_abort(
      "All elements of {.arg at} must be named (use {.code list(factor = level, cov = value)})."
    )
  }
  cov_names <- fit_obj$param_info$continuous_covariates %||% character(0)
  all_factors <- fit_obj$param_info$factors %||% character(0)
  all_factors <- all_factors[nzchar(all_factors) & !is.na(all_factors)]
  valid_names <- c(all_factors, cov_names)
  bad_names <- setdiff(names(at), valid_names)
  if (length(bad_names) > 0L) {
    cli::cli_abort(c(
      "Unknown name{?s} in {.arg at}: {.field {bad_names}}.",
      "i" = "Valid names are the fit's factors and continuous covariates: {.field {valid_names}}.",
      "x" = "Did you mistype a factor or covariate name?"
    ))
  }
  data_used <- fit_obj$data
  for (nm in names(at)) {
    v <- at[[nm]]
    if (length(v) < 1L) {
      cli::cli_abort(c(
        "{.field {nm}} has length 0.",
        "i" = "Each {.arg at} entry must be a non-empty vector."
      ))
    }
    if (nm %in% all_factors) {
      observed <- sort(unique(as.character(data_used[[nm]])))
      bad_vals <- setdiff(as.character(v), observed)
      if (length(bad_vals) > 0L) {
        cli::cli_abort(c(
          "{.field {nm}} = {.val {bad_vals}} not an observed level.",
          "i" = "Observed levels: {.val {observed}}."
        ))
      }
    } else {
      v_num <- suppressWarnings(as.numeric(v))
      if (any(is.na(v_num)) || any(!is.finite(v_num))) {
        cli::cli_abort(c(
          "{.field {nm}} value{?s} {.val {as.character(v)}} not finite numeric.",
          "i" = "Continuous-covariate {.arg at} entries must be a single finite numeric value."
        ))
      }
      if (length(v) > 1L) {
        cli::cli_warn(c(
          "{.arg at${nm}} has length {length(v)}; using first value {.val {v_num[1]}}.",
          "i" = "Pass a single numeric value per continuous covariate."
        ))
      }
    }
  }
  invisible(NULL)
}

# Build the conditioned reference grid (level_combos) and the averaging-matrix
# design (ref_X = A %*% X_full) for k EMMs/comparisons. Single-parameter port of
# beezdemand .tmb_build_emm_ref_grid() (R/tmb-methods.R:3186); the only fitted
# design is fit_obj$formula_details$X and the only factor set is
# fit_obj$param_info$factors.
.dd_build_emm_ref_grid <- function(
  fit_obj,
  at = NULL,
  factors_in_emm = NULL,
  validate = TRUE
) {
  cov_names <- fit_obj$param_info$continuous_covariates %||% character(0)
  fitted_factors <- fit_obj$param_info$factors %||% character(0)
  fitted_factors <- fitted_factors[nzchar(fitted_factors) & !is.na(fitted_factors)]

  use_factors <- fitted_factors
  if (!is.null(factors_in_emm)) {
    if (length(factors_in_emm) == 0L) {
      use_factors <- character(0)            # ~ 1: marginalize everything
    } else {
      use_factors <- .dd_resolve_retained_factors(factors_in_emm, fitted_factors)
    }
  }

  if (isTRUE(validate)) .dd_validate_at(fit_obj, at)

  retained_factors <- use_factors
  is_intercept_only <- length(fitted_factors) == 0L && length(cov_names) == 0L

  if (is_intercept_only) {
    return(list(
      level_combos = NULL, ref_X = NULL, use_factors = character(0),
      cov_names = character(0), is_intercept_only = TRUE
    ))
  }

  data_used <- fit_obj$data

  factor_level_set <- function(f) {
    lv <- levels(data_used[[f]])
    if (is.null(lv)) lv <- sort(unique(as.character(data_used[[f]])))
    if (!is.null(at) && f %in% names(at)) lv <- lv[lv %in% as.character(at[[f]])]
    lv
  }
  fitted_levels <- stats::setNames(
    lapply(fitted_factors, factor_level_set), fitted_factors
  )
  if (length(fitted_factors) > 0L &&
      any(vapply(fitted_levels, length, integer(1)) == 0L)) {
    cli::cli_abort(c(
      "{.arg at} filter produced an empty reference grid.",
      "i" = "Check that the supplied factor levels exist in the data."
    ))
  }

  make_key <- function(df, cols) {
    if (length(cols) == 0L) return(rep("", nrow(df)))
    do.call(paste, c(lapply(cols, function(cc) as.character(df[[cc]])),
                     list(sep = "\r")))
  }
  as_training_factor <- function(values, f) {
    factor(values, levels = levels(data_used[[f]]) %||%
             sort(unique(as.character(data_used[[f]]))))
  }

  # Full factorial grid over ALL fitted factors (equal-weight averaging target).
  if (length(fitted_factors) > 0L) {
    full_combos <- do.call(expand.grid, c(
      lapply(fitted_factors, function(f) as_training_factor(fitted_levels[[f]], f)),
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
    names(full_combos) <- fitted_factors
  } else {
    full_combos <- data_used[1L, integer(0), drop = FALSE]
  }

  # Retained reference grid: crossing of retained factors, ordered by level
  # index, filtered to OBSERVED combinations (semi_join analog).
  if (length(retained_factors) > 0L) {
    level_combos <- do.call(expand.grid, c(
      lapply(retained_factors, function(f) as_training_factor(fitted_levels[[f]], f)),
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
    names(level_combos) <- retained_factors
    ord <- do.call(order, lapply(retained_factors,
                                 function(f) as.integer(level_combos[[f]])))
    level_combos <- level_combos[ord, , drop = FALSE]
    observed_keys <- make_key(
      unique(data_used[, retained_factors, drop = FALSE]), retained_factors
    )
    level_combos <- level_combos[
      make_key(level_combos, retained_factors) %in% observed_keys, , drop = FALSE
    ]
    if (nrow(level_combos) == 0L) {
      cli::cli_abort(c(
        "{.arg at} filter produced an empty reference grid.",
        "i" = "Check that the supplied factor levels are not mutually exclusive."
      ))
    }
    rownames(level_combos) <- NULL
  } else {
    level_combos <- data_used[1L, integer(0), drop = FALSE]
  }

  # Continuous covariates: training mean unless overridden via `at`.
  if (length(cov_names) > 0L) {
    for (cv in cov_names) {
      cv_value <- mean(data_used[[cv]], na.rm = TRUE)
      if (!is.null(at) && cv %in% names(at)) cv_value <- as.numeric(at[[cv]][1])
      full_combos[[cv]] <- cv_value
      level_combos[[cv]] <- cv_value
    }
  }

  # Pin the rebuilt basis to the FITTED design's contrasts, then verify (and
  # reorder to) the fitted column set — abort loudly on mismatch.
  fitted_X <- fit_obj$formula_details$X
  X_full <- stats::model.matrix(
    stats::as.formula(build_fixed_rhs(
      factors = fitted_factors,
      factor_interaction = fit_obj$param_info$factor_interaction,
      continuous_covariates = cov_names,
      data = data_used
    )),
    data = full_combos,
    contrasts.arg = attr(fitted_X, "contrasts")
  )
  fitted_cols <- colnames(fitted_X)
  if (!is.null(fitted_cols)) {
    if (!setequal(colnames(X_full), fitted_cols)) {
      cli::cli_abort(c(
        "Could not reproduce the fitted design matrix for the EMM grid.",
        "i" = "Rebuilt columns: {.val {colnames(X_full)}}.",
        "i" = "Fitted columns: {.val {fitted_cols}}.",
        "x" = "This can happen if the model's factor levels or contrasts changed after fitting."
      ))
    }
    X_full <- X_full[, fitted_cols, drop = FALSE]
  }

  # Averaging matrix A (n_retained x n_full): equal weight 1/m on the m
  # full-grid rows matching each retained cell. ref_X = A %*% X_full keeps
  # ncol == length(beta_k).
  full_keys <- make_key(full_combos, retained_factors)
  ret_keys <- make_key(level_combos, retained_factors)
  A <- matrix(0, nrow = nrow(level_combos), ncol = nrow(full_combos))
  for (r in seq_len(nrow(level_combos))) {
    sel <- which(full_keys == ret_keys[r])
    A[r, sel] <- 1 / length(sel)
  }
  ref_X <- A %*% X_full
  colnames(ref_X) <- colnames(X_full)

  list(
    level_combos = level_combos,
    ref_X = ref_X,
    use_factors = retained_factors,
    cov_names = cov_names,
    is_intercept_only = FALSE
  )
}
```

  Also add the null-coalescing helper if not already defined in the package:

```r
# R/dd-comparisons.R (top, if `%||%` is not already exported/imported elsewhere)
`%||%` <- function(x, y) if (is.null(x)) y else x
```

- [ ] **Step 4: Run the test, expect pass.**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: the three `.dd_build_emm_ref_grid()` `it()` blocks pass.

- [ ] **Step 5: Commit.**

  Run:
  ```
  git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-comparisons.R tests/testthat/test-dd-comparisons.R
  git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(emmeans): add .dd_build_emm_ref_grid averaging-matrix helper

Single-parameter port of beezdemand .tmb_build_emm_ref_grid for the
discounting log-k design: full-factorial grid, fitted-contrasts basis,
averaging matrix A, ref_X = A %*% X_full.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task E.2: `get_dd_param_emms()` — marginal means for k

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-comparisons.R`
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/NAMESPACE` (via roxygen `@export`)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R`

- [ ] **Step 1: Write the failing test (non-circular: recompute marginal means from coefficients %*% model.matrix).**

```r
# tests/testthat/test-dd-comparisons.R  (append)
describe("get_dd_param_emms()", {

  it("returns the contract columns and back-transforms log k to k = exp(k_log)", {
    set.seed(10)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.8, seed = 10
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm <- get_dd_param_emms(fit)
    expect_named(emm, c("level", "k", "k_log", "std.error", "conf.low", "conf.high"))
    expect_equal(nrow(emm), nlevels(fit$data$condition))
    expect_equal(emm$k, exp(emm$k_log), tolerance = 1e-12)
    # CIs back-transform from the log scale (so always positive, bracket k).
    expect_true(all(emm$conf.low > 0))
    expect_true(all(emm$conf.low <= emm$k & emm$k <= emm$conf.high))
  })

  it("recomputes per-cell log-k means NON-CIRCULARLY from coef %*% model.matrix", {
    set.seed(11)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.8, seed = 11
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm <- get_dd_param_emms(fit)

    # Truth, by hand: beta_k %*% (per-cell design row), never via emmeans.
    beta <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"])
    lv <- levels(fit$data$condition)
    nd <- data.frame(condition = factor(lv, levels = lv))
    Xc <- stats::model.matrix(~ condition, data = nd,
                              contrasts.arg = attr(fit$formula_details$X, "contrasts"))
    Xc <- Xc[, colnames(fit$formula_details$X), drop = FALSE]
    expected_log_k <- as.numeric(Xc %*% beta)
    expect_equal(emm$k_log, expected_log_k, tolerance = 1e-8)
  })

  it("intercept-only fit returns a single (Intercept) row equal to exp(beta_k[1])", {
    set.seed(12)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 40, family = "sltb", equation = "mazur", seed = 12
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    emm <- get_dd_param_emms(fit)
    expect_equal(nrow(emm), 1L)
    expect_equal(emm$level, "(Intercept)")
    beta1 <- unname(fit$model$coefficients[names(fit$model$coefficients) == "beta_k"][1])
    expect_equal(emm$k_log, beta1, tolerance = 1e-10)
    expect_equal(emm$k, exp(beta1), tolerance = 1e-10)
  })

  it("marginalizes an omitted factor: factors_in_emm = character(0) gives one cell", {
    set.seed(13)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.8, seed = 13
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm0 <- get_dd_param_emms(fit, factors_in_emm = character(0))
    expect_equal(nrow(emm0), 1L)
    # The grand mean of the per-cell log-k EMMs (equal weights).
    emm <- get_dd_param_emms(fit)
    expect_equal(emm0$k_log, mean(emm$k_log), tolerance = 1e-8)
  })
})
```

- [ ] **Step 2: Run the test, expect failure (`get_dd_param_emms` undefined).**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: `could not find function "get_dd_param_emms"`.

- [ ] **Step 3: Implement `get_dd_param_emms()` in `R/dd-comparisons.R`.**

  Simplification of beezdemand `get_demand_param_emms.beezdemand_tmb()`
  (`R/tmb-methods.R:3445-3579`): the only beta block is `beta_k`; the only
  design matrix is `fit$formula_details$X`; there is no `param` switch and no
  generic/method split (discounting exposes a plain function, not an S3
  generic). Column names follow the PINNED CONTRACT
  (`level, k, k_log, std.error, conf.low, conf.high`) — note `estimate`→`k`,
  `estimate_log`→`k_log` versus beezdemand.

```r
# R/dd-comparisons.R  (append)

#' Estimated marginal means of the discount rate \code{k}
#'
#' @description
#' Computes estimated marginal means (EMMs) of the discount rate \code{k} from a
#' fitted \code{beezdiscounting_tmb} model. EMMs are computed on the
#' \code{log k} scale (linear in the fixed-effect coefficients) using the
#' averaging-matrix reference grid, then back-transformed with \code{exp()} so
#' that \code{k = exp(k_log)}. Standard errors use the \code{beta_k} block of
#' \code{TMB::sdreport()}'s fixed-effect covariance; intervals are Wald on the
#' log scale and exponentiated.
#'
#' @param fit A \code{beezdiscounting_tmb} object.
#' @param factors_in_emm Character vector of factors to retain in the EMM
#'   reference grid. A strict subset marginalizes the omitted factors with equal
#'   weights across the full crossing of their levels (emmeans' default
#'   \code{weights = "equal"}); \code{NULL} (default) retains all fitted factors;
#'   \code{character(0)} marginalizes everything to a single grand-mean cell.
#' @param at Named list specifying factor levels and/or continuous-covariate
#'   values for conditional EMMs (one numeric per covariate; multiple values
#'   warn and use the first). \code{at} on an omitted factor restricts the
#'   level set averaged over.
#' @param ci_level Numeric confidence level for intervals (default 0.95).
#' @param ... Additional arguments (currently unused).
#'
#' @return A tibble with columns \code{level}, \code{k}, \code{k_log},
#'   \code{std.error}, \code{conf.low}, \code{conf.high}. \code{k_log} is the
#'   marginal mean on the \code{log k} scale; \code{k = exp(k_log)};
#'   \code{std.error} is the SE of \code{k_log}; the intervals are on the
#'   \code{k} (natural) scale.
#'
#' @export
get_dd_param_emms <- function(
  fit,
  factors_in_emm = NULL,
  at = NULL,
  ci_level = 0.95,
  ...
) {
  coefs <- fit$model$coefficients
  beta_idx <- which(names(coefs) == "beta_k")
  beta <- unname(coefs[beta_idx])

  # beta_k block of the fixed-effect covariance (sdr$cov.fixed), aligned by the
  # opt$par naming. Falls back to a diagonal of squared SEs if sdreport failed.
  vcov_mat <- NULL
  sdr <- fit$sdr
  if (!is.null(sdr) && !is.null(sdr$cov.fixed)) {
    full_vcov <- as.matrix(sdr$cov.fixed)
    par_names <- names(fit$opt$par)
    target_idx <- which(par_names == "beta_k")
    if (length(target_idx) == length(beta)) {
      vcov_mat <- full_vcov[target_idx, target_idx, drop = FALSE]
    }
  }
  if (is.null(vcov_mat)) {
    se_vals <- fit$model$se[beta_idx]
    vcov_mat <- diag(se_vals^2, nrow = length(se_vals))
  }

  .dd_validate_at(fit, at)

  grid <- .dd_build_emm_ref_grid(
    fit, at = at, factors_in_emm = factors_in_emm, validate = FALSE
  )
  z <- stats::qnorm((1 + ci_level) / 2)

  if (isTRUE(grid$is_intercept_only)) {
    est <- beta[1L]
    se <- sqrt(vcov_mat[1L, 1L])
    return(tibble::tibble(
      level = "(Intercept)",
      k = exp(est),
      k_log = est,
      std.error = se,
      conf.low = exp(est - z * se),
      conf.high = exp(est + z * se)
    ))
  }

  use_factors <- grid$use_factors
  cov_names <- grid$cov_names
  level_combos <- grid$level_combos
  ref_X <- grid$ref_X

  if (ncol(ref_X) != length(beta)) {
    cli::cli_abort(c(
      "Reference-grid design has {ncol(ref_X)} column{?s} but the fitted \\
       coefficient vector has {length(beta)}.",
      "x" = "Design basis mismatch; cannot evaluate EMMs."
    ))
  }

  level_label_for <- function(i) {
    if (length(use_factors) > 0L) {
      paste(vapply(use_factors, function(f)
        paste0(f, "=", level_combos[[f]][i]), character(1)), collapse = ", ")
    } else if (length(cov_names) > 0L) {
      paste(vapply(cov_names, function(cv)
        paste0(cv, "=", level_combos[[cv]][i]), character(1)), collapse = ", ")
    } else {
      "(Intercept)"
    }
  }

  cell_est <- as.numeric(ref_X %*% beta)
  cell_se <- sqrt(diag(ref_X %*% vcov_mat %*% t(ref_X)))

  tibble::tibble(
    level = vapply(seq_len(nrow(ref_X)), level_label_for, character(1)),
    k = exp(cell_est),
    k_log = cell_est,
    std.error = cell_se,
    conf.low = exp(cell_est - z * cell_se),
    conf.high = exp(cell_est + z * cell_se)
  )
}
```

- [ ] **Step 4: Run the test, expect pass.**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: all `get_dd_param_emms()` `it()` blocks pass.

- [ ] **Step 5: Regenerate NAMESPACE and commit.**

  Run:
  ```
  Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting")'
  git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-comparisons.R tests/testthat/test-dd-comparisons.R NAMESPACE man/get_dd_param_emms.Rd
  git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(emmeans): add get_dd_param_emms for the discount rate k

Marginal means of k on the log scale via the averaging-matrix grid, SE
from the beta_k block of sdr\$cov.fixed, exp() back-transform. Returns the
contract columns level, k, k_log, std.error, conf.low, conf.high.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task E.3: `get_dd_comparisons()` — contrasts on log k (log10 + ratios)

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-comparisons.R`
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/NAMESPACE` (via roxygen `@export`)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R`

- [ ] **Step 1: Write the failing test (EMM↔contrast log(10) invariance to 1e-8; choose(n,2) rows; adjust validation).**

```r
# tests/testthat/test-dd-comparisons.R  (append)
describe("get_dd_comparisons()", {

  it("rejects emmeans-only p-adjust methods (must be in stats::p.adjust.methods)", {
    set.seed(20)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.8, seed = 20
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    expect_error(get_dd_comparisons(fit, adjust = "tukey"), "p-value adjustment")
    expect_error(get_dd_comparisons(fit, adjust = "sidak"), "p-value adjustment")
  })

  it("returns a beezdiscounting_comparison with k$emmeans/contrasts_log10/contrasts_ratio", {
    set.seed(21)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 60, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.8, seed = 21
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit)
    expect_s3_class(res, "beezdiscounting_comparison")
    expect_named(res, "k")
    expect_named(res$k, c("emmeans", "contrasts_log10", "contrasts_ratio"))
    # 2 levels -> choose(2,2) = 1 pairwise contrast.
    expect_equal(nrow(res$k$contrasts_log10), 1L)
    expect_equal(nrow(res$k$contrasts_ratio), 1L)
  })

  it("produces choose(n,2) pairwise rows for an n-level factor", {
    set.seed(22)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 90, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = 0.6, seed = 22
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit)
    expect_equal(nrow(res$k$contrasts_log10), choose(3L, 2L))
  })

  it("trt.vs.ctrl gives n-1 rows (each level vs the reference)", {
    set.seed(23)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 90, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = 0.6, seed = 23
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit, contrast_type = "trt.vs.ctrl")
    expect_equal(nrow(res$k$contrasts_log10), 3L - 1L)
  })

  it("EMM<->contrast invariance: log10 estimate == (k_log[i]-k_log[j])/log(10) to 1e-8", {
    set.seed(24)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 90, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = 0.6, seed = 24
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm <- get_dd_param_emms(fit)
    res <- get_dd_comparisons(fit)

    # The contrasts are emmeans rows ordered by combn(seq_len(n), 2): (1,2),(1,3),(2,3).
    cmb <- utils::combn(nrow(emm), 2L)
    expected_log10 <- (emm$k_log[cmb[1L, ]] - emm$k_log[cmb[2L, ]]) / log(10)
    expect_equal(res$k$contrasts_log10$estimate, expected_log10, tolerance = 1e-8)

    # Ratio block: exp(est_log) == 10^(log10 estimate).
    expect_equal(res$k$contrasts_ratio$ratio,
                 10^res$k$contrasts_log10$estimate, tolerance = 1e-8)
  })

  it("contrast_by adjusts p-values per by-cell", {
    set.seed(25)
    # Two crossed factors so we can condition contrasts of one by the other.
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 120, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.7, seed = 25
    )
    dat$grp <- factor(rep(c("A", "B"), length.out = nrow(dat)))
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = c("condition", "grp"), factor_interaction = TRUE,
                      multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit, compare_specs = ~ condition + grp,
                              contrast_by = "grp", adjust = "holm")
    # One condition contrast within each of the 2 grp by-cells.
    expect_true(nrow(res$k$contrasts_log10) == 2L)
    expect_true("grp" %in% names(res$k$contrasts_log10))
  })
})
```

- [ ] **Step 2: Run the test, expect failure (`get_dd_comparisons` undefined).**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: `could not find function "get_dd_comparisons"`.

- [ ] **Step 3: Implement `get_dd_comparisons()` + `.dd_compare_block()` in `R/dd-comparisons.R`.**

  Single-parameter port of beezdemand `get_demand_comparisons.beezdemand_tmb()`
  + `.tmb_compare_one_param()` (`R/tmb-methods.R:3653-4074`). Adaptations:
  remove the `param` argument and the per-parameter `lapply` (one parameter,
  `k`); the result list is keyed by the single name `"k"`; drop all
  collapse-aliasing / per-param `by_map` machinery (no `collapse_levels` in the
  discounting model — `contrast_by` resolves by direct name match only); class
  is `beezdiscounting_comparison`. The contrast math is verbatim:
  `dx = ref_X[i,] - ref_X[j,]`, `est_log = sum(dx*beta)`,
  `se_log = sqrt(t(dx) %*% vcov %*% dx)`, `z = est_log/se_log`,
  `p.adjust(2*pnorm(-abs(z)))` per block, `contrasts_log10 = est_log/log(10)`,
  `contrasts_ratio = exp(est_log)` with `exp(±z*se_log)` intervals.

```r
# R/dd-comparisons.R  (append)

#' Pairwise / treatment-vs-control contrasts of the discount rate \code{k}
#'
#' @description
#' Computes factor-level contrasts of \code{k} from a \code{beezdiscounting_tmb}
#' model. Contrasts are formed on the \code{log k} scale via the
#' averaging-matrix reference grid and reported two ways: \code{contrasts_log10}
#' (differences on the \code{log10 k} scale, with asymptotic-\emph{z} Wald
#' intervals and p-values adjusted by \code{adjust}) and, when
#' \code{report_ratios = TRUE}, \code{contrasts_ratio} (multiplicative
#' \code{exp(log-difference)} ratios).
#'
#' @param fit A \code{beezdiscounting_tmb} object.
#' @param compare_specs Optional one-sided formula naming the factor subset to
#'   contrast (e.g. \code{~ condition}). Omitted fitted factors are marginalized
#'   with equal weights. \code{NULL} (default) retains all fitted factors.
#' @param contrast_type \code{"pairwise"} (all pairs, factor-level order) or
#'   \code{"trt.vs.ctrl"} (each level vs. the reference/first level).
#' @param contrast_by Optional \code{NULL} (default) or character vector of
#'   factor name(s) within \code{compare_specs} to condition contrasts on.
#'   Within each observed by-cell, contrasts are computed over the remaining
#'   factors with p-value adjustment applied \strong{per by-cell}. A
#'   \code{contrast_by} factor absent from \code{compare_specs} aborts.
#' @param adjust P-value adjustment method; must be one of
#'   \code{stats::p.adjust.methods} (default \code{"holm"}). emmeans-only
#'   methods (\code{"tukey"}, \code{"sidak"}, ...) are rejected.
#' @param at Named list conditioning factor levels / covariate values, as in
#'   \code{\link{get_dd_param_emms}}.
#' @param ci_level Numeric confidence level (default 0.95).
#' @param report_ratios Logical; include the \code{contrasts_ratio} block
#'   (default \code{TRUE}).
#' @param ... Additional arguments (reserved).
#'
#' @return An object of class \code{"beezdiscounting_comparison"}: a list with a
#'   single element \code{k}, itself a list of \code{emmeans} (cell means),
#'   \code{contrasts_log10} (\code{contrast}, \code{estimate}, \code{std.error},
#'   \code{statistic}, \code{df}, \code{conf.low}, \code{conf.high},
#'   \code{p.value}; on the \code{log10 k} scale), and (if \code{report_ratios})
#'   \code{contrasts_ratio}. When \code{contrast_by} is active the contrast
#'   tables gain leading by-column(s). Attributes \code{backend},
#'   \code{adjustment_method}, \code{compare_specs_used}, \code{contrast_type_used},
#'   and \code{contrast_by_used} describe the call.
#'
#' @export
get_dd_comparisons <- function(
  fit,
  compare_specs = NULL,
  contrast_type = c("pairwise", "trt.vs.ctrl"),
  contrast_by = NULL,
  adjust = "holm",
  at = NULL,
  ci_level = 0.95,
  report_ratios = TRUE,
  ...
) {
  contrast_type <- match.arg(contrast_type)

  if (!isTRUE(adjust %in% stats::p.adjust.methods)) {
    cli::cli_abort(c(
      "{.arg adjust} = {.val {adjust}} is not a valid p-value adjustment method.",
      "i" = "Valid methods: {.val {stats::p.adjust.methods}}.",
      "x" = "emmeans-only methods (e.g. {.val tukey}, {.val sidak}) are unavailable here (asymptotic z + {.fn stats::p.adjust})."
    ))
  }

  fitted_factors <- fit$param_info$factors %||% character(0)
  fitted_factors <- fitted_factors[nzchar(fitted_factors) & !is.na(fitted_factors)]

  # contrast_by boundary validation (typos error here, once).
  if (!is.null(contrast_by)) {
    if (!is.character(contrast_by)) {
      cli::cli_abort("{.arg contrast_by} must be {.code NULL} or a character vector of factor name(s).")
    }
    if (length(contrast_by) == 0L) {
      contrast_by <- NULL
    } else {
      bad_by <- setdiff(contrast_by, fitted_factors)
      if (length(bad_by) > 0L) {
        cli::cli_abort(c(
          "{.arg contrast_by} names factor{?s} not in the fit: {.val {bad_by}}.",
          "i" = "Fitted factors: {.val {fitted_factors}}."
        ))
      }
    }
  }

  # Resolve the retained factor set from compare_specs (canonical).
  factors_in_emm <- NULL
  if (!is.null(compare_specs)) {
    if (!inherits(compare_specs, "formula")) {
      cli::cli_abort("{.arg compare_specs} must be a one-sided formula (e.g. {.code ~ condition}).")
    }
    factors_in_emm <- all.vars(compare_specs)
    bad <- setdiff(factors_in_emm, fitted_factors)
    if (length(bad) > 0L) {
      cli::cli_abort(c(
        "{.arg compare_specs} names factor{?s} not in the fit: {.val {bad}}.",
        "i" = "Fitted factors: {.val {fitted_factors}}."
      ))
    }
  }

  .dd_validate_at(fit, at)

  block <- .dd_compare_block(
    fit, factors_in_emm, contrast_type, adjust, at, ci_level,
    report_ratios, contrast_by
  )

  results_list <- list(k = block)
  class(results_list) <- "beezdiscounting_comparison"
  attr(results_list, "backend") <- "tmb"
  attr(results_list, "compare_specs_used") <- if (is.null(compare_specs)) {
    "all fitted factors"
  } else {
    deparse(compare_specs)
  }
  attr(results_list, "contrast_type_used") <- contrast_type
  any_by_applied <- !is.null(attr(block, "effective_by")) &&
    length(attr(block, "effective_by")) > 0L
  attr(results_list, "contrast_by_used") <-
    if (is.null(contrast_by) || !any_by_applied) "NULL" else paste(contrast_by, collapse = ", ")
  attr(results_list, "adjustment_method") <- adjust
  results_list
}

# Build the single-parameter (k) comparison block. Returns
# list(emmeans, contrasts_log10[, contrasts_ratio]); contrasts_log10 carries an
# `std_labels` attribute (structured ref-grid labels) for the tidy method, and
# the block carries an `effective_by` attribute for the public-function metadata.
.dd_compare_block <- function(fit, factors_in_emm, contrast_type, adjust,
                              at, ci_level, report_ratios, contrast_by = NULL) {
  coefs <- fit$model$coefficients
  beta <- unname(coefs[names(coefs) == "beta_k"])

  vcov_mat <- NULL
  sdr <- fit$sdr
  if (!is.null(sdr) && !is.null(sdr$cov.fixed)) {
    full_vcov <- as.matrix(sdr$cov.fixed)
    par_names <- names(fit$opt$par)
    target_idx <- which(par_names == "beta_k")
    if (length(target_idx) == length(beta)) {
      vcov_mat <- full_vcov[target_idx, target_idx, drop = FALSE]
    }
  }
  if (is.null(vcov_mat)) {
    se_vals <- fit$model$se[names(coefs) == "beta_k"]
    vcov_mat <- diag(se_vals^2, nrow = length(se_vals))
  }

  grid <- .dd_build_emm_ref_grid(
    fit, at = at, factors_in_emm = factors_in_emm, validate = FALSE
  )
  z <- stats::qnorm((1 + ci_level) / 2)
  ln10 <- log(10)

  empty_log10 <- tibble::tibble(
    contrast = character(), estimate = numeric(), std.error = numeric(),
    statistic = numeric(), df = numeric(),
    conf.low = numeric(), conf.high = numeric(), p.value = numeric()
  )
  empty_ratio <- tibble::tibble(
    contrast = character(), ratio = numeric(),
    conf.low = numeric(), conf.high = numeric(), p.value = numeric()
  )
  finish_empty <- function(emm_block) {
    out <- list(emmeans = emm_block, contrasts_log10 = empty_log10)
    attr(out$contrasts_log10, "std_labels") <- character()
    if (report_ratios) out$contrasts_ratio <- empty_ratio
    attr(out, "effective_by") <- character(0)
    out
  }

  if (isTRUE(grid$is_intercept_only)) {
    est <- beta[1L]
    se <- sqrt(vcov_mat[1L, 1L])
    emm_block <- tibble::tibble(
      level = "(Intercept)", k = exp(est), k_log = est, std.error = se,
      conf.low = exp(est - z * se), conf.high = exp(est + z * se)
    )
    return(finish_empty(emm_block))
  }

  use_factors <- grid$use_factors
  cov_names <- grid$cov_names
  level_combos <- grid$level_combos
  ref_X <- grid$ref_X
  n <- nrow(ref_X)

  cell_est <- as.numeric(ref_X %*% beta)
  cell_se <- sqrt(diag(ref_X %*% vcov_mat %*% t(ref_X)))

  native_label_f <- function(i, fs) {
    if (length(fs) > 0L) {
      paste(vapply(fs, function(f)
        paste0(f, "=", as.character(level_combos[[f]][i])), character(1)),
        collapse = ", ")
    } else if (length(cov_names) > 0L) {
      paste(vapply(cov_names, function(cv)
        paste0(cv, "=", level_combos[[cv]][i]), character(1)), collapse = ", ")
    } else {
      "(Intercept)"
    }
  }
  std_label_f <- function(i, fs) {
    if (length(fs) > 0L) {
      paste(vapply(fs, function(f)
        as.character(level_combos[[f]][i]), character(1)), collapse = " ")
    } else {
      native_label_f(i, fs)
    }
  }
  native_label <- function(i) native_label_f(i, use_factors)

  emm_block <- tibble::tibble(
    level = vapply(seq_len(n), native_label, character(1)),
    k = exp(cell_est), k_log = cell_est, std.error = cell_se,
    conf.low = exp(cell_est - z * cell_se),
    conf.high = exp(cell_est + z * cell_se)
  )

  # Resolve contrast_by for the (single) k parameter: direct name match only
  # (no collapse aliasing in the discounting model). Each by-var must be in the
  # retained set; abort otherwise.
  effective_by <- character(0)
  if (n >= 2L && !is.null(contrast_by)) {
    effective_by <- intersect(contrast_by, use_factors)
    not_in <- setdiff(contrast_by, use_factors)
    if (length(not_in) > 0L) {
      cli::cli_abort(c(
        "{cli::qty(not_in)}{.arg contrast_by} factor{?s} {.val {not_in}} {?is/are} not in {.arg compare_specs}.",
        "i" = "{cli::qty(use_factors)}{.arg compare_specs} factor{?s}: {.val {use_factors}}.",
        "x" = "Name the by-variable(s) in {.arg compare_specs} to condition contrasts on them."
      ))
    }
    # Redundant-by: a single retained factor equal to the by-set means there is
    # nothing left to contrast within a cell; ignore by-grouping.
    if (length(use_factors) == 1L && identical(sort(use_factors), sort(effective_by))) {
      message("  `contrast_by` is redundant with `compare_specs`; ignoring it.")
      effective_by <- character(0)
    }
  }

  comparison_factors <- setdiff(use_factors, effective_by)

  if (length(effective_by) > 0L) {
    by_key <- do.call(paste, c(
      lapply(effective_by, function(f) as.character(level_combos[[f]])),
      list(sep = "\r")
    ))
    blocks <- lapply(unique(by_key), function(k) which(by_key == k))
  } else {
    blocks <- list(seq_len(n))
  }

  do_block <- function(rows) {
    m <- length(rows)
    if (m < 2L) return(NULL)
    if (contrast_type == "pairwise") {
      cmb <- utils::combn(m, 2L)
      lhs <- rows[cmb[1L, ]]
      rhs <- rows[cmb[2L, ]]
    } else {
      lhs <- rows[seq.int(2L, m)]
      rhs <- rep(rows[1L], m - 1L)
    }
    est_log <- numeric(length(lhs))
    se_log <- numeric(length(lhs))
    native <- character(length(lhs))
    stdlab <- character(length(lhs))
    for (k in seq_along(lhs)) {
      dx <- ref_X[lhs[k], ] - ref_X[rhs[k], ]
      est_log[k] <- sum(dx * beta)
      se_log[k] <- sqrt(as.numeric(t(dx) %*% vcov_mat %*% dx))
      native[k] <- paste(native_label_f(lhs[k], comparison_factors), "-",
                         native_label_f(rhs[k], comparison_factors))
      stdlab[k] <- paste(std_label_f(lhs[k], comparison_factors), "-",
                         std_label_f(rhs[k], comparison_factors))
    }
    zstat <- est_log / se_log
    p_adj <- stats::p.adjust(2 * stats::pnorm(-abs(zstat)), method = adjust)
    est_log10 <- est_log / ln10
    se_log10 <- se_log / ln10
    list(
      log10 = tibble::tibble(
        contrast = native, estimate = est_log10, std.error = se_log10,
        statistic = zstat, df = Inf,
        conf.low = est_log10 - z * se_log10,
        conf.high = est_log10 + z * se_log10, p.value = p_adj
      ),
      ratio = tibble::tibble(
        contrast = native, ratio = exp(est_log),
        conf.low = exp(est_log - z * se_log),
        conf.high = exp(est_log + z * se_log), p.value = p_adj
      ),
      std_labels = stdlab,
      first_row = rows[1L]
    )
  }

  by_cols_for <- function(first_row, nrows) {
    if (length(effective_by) == 0L) return(NULL)
    cols <- lapply(effective_by, function(f) {
      rep(as.character(level_combos[[f]][first_row]), nrows)
    })
    tibble::as_tibble(stats::setNames(cols, effective_by))
  }

  block_results <- Filter(Negate(is.null), lapply(blocks, do_block))
  if (length(block_results) == 0L) {
    res <- finish_empty(emm_block)
    attr(res, "effective_by") <- effective_by
    return(res)
  }

  log10_parts <- lapply(block_results, function(r) {
    bc <- by_cols_for(r$first_row, nrow(r$log10))
    if (is.null(bc)) r$log10 else dplyr::bind_cols(bc, r$log10)
  })
  ratio_parts <- lapply(block_results, function(r) {
    bc <- by_cols_for(r$first_row, nrow(r$ratio))
    if (is.null(bc)) r$ratio else dplyr::bind_cols(bc, r$ratio)
  })

  contrasts_log10 <- dplyr::bind_rows(log10_parts)
  attr(contrasts_log10, "std_labels") <- unlist(
    lapply(block_results, function(r) r$std_labels), use.names = FALSE
  )

  out <- list(emmeans = emm_block, contrasts_log10 = contrasts_log10)
  if (report_ratios) {
    contrasts_ratio <- dplyr::bind_rows(ratio_parts)
    attr(contrasts_ratio, "std_labels") <- attr(contrasts_log10, "std_labels")
    out$contrasts_ratio <- contrasts_ratio
  }
  attr(out, "effective_by") <- effective_by
  out
}
```

- [ ] **Step 4: Run the test, expect pass.**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: all `get_dd_comparisons()` `it()` blocks pass (invariance to 1e-8, choose(n,2) rows, trt.vs.ctrl n-1 rows, adjust rejection, contrast_by per-cell).

- [ ] **Step 5: Regenerate NAMESPACE and commit.**

  Run:
  ```
  Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting")'
  git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-comparisons.R tests/testthat/test-dd-comparisons.R NAMESPACE man/get_dd_comparisons.Rd
  git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(emmeans): add get_dd_comparisons for k contrasts (log10 + ratios)

Pairwise / trt.vs.ctrl contrasts of log k via dx = ref_X[i,]-ref_X[j,],
asymptotic-z Wald, stats::p.adjust per contrast_by cell, log10 + ratio
blocks; classed beezdiscounting_comparison.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```

---

### Task E.4: `tidy.beezdiscounting_comparison()` — flat cross-backend frame

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-comparisons.R`
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/NAMESPACE` (via roxygen `@export` + `@importFrom generics tidy`)
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/DESCRIPTION` (`Imports: generics`)
- Test: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R`

- [ ] **Step 1: Write the failing test (flat columns; structured std labels; exponentiate ratios; by-columns).**

```r
# tests/testthat/test-dd-comparisons.R  (append)
describe("tidy.beezdiscounting_comparison()", {

  it("flattens to the cross-backend column set (param='k')", {
    set.seed(30)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 90, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = 0.6, seed = 30
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit)
    td <- generics::tidy(res)
    expect_equal(names(td),
                 c("param", "contrast", "estimate", "std.error", "statistic",
                   "df", "conf.low", "conf.high", "p.value"))
    expect_true(all(td$param == "k"))
    expect_equal(nrow(td), choose(3L, 2L))
    # statistic is asymptotic z with df = Inf on the TMB backend.
    expect_true(all(is.infinite(td$df)))
  })

  it("uses structured std_labels for the contrast column (level values, not k=..)", {
    set.seed(31)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 90, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = 0.6, seed = 31
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit)
    td <- generics::tidy(res)
    # labels are "lvl_i - lvl_j", never the native "condition=.. - condition=.."
    expect_false(any(grepl("condition=", td$contrast)))
    expect_true(all(grepl(" - ", td$contrast)))
  })

  it("EMM<->tidy invariance holds end to end (td$estimate == (k_log[i]-k_log[j])/log(10))", {
    set.seed(32)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 90, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = 0.6, seed = 32
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    emm <- get_dd_param_emms(fit)
    td <- generics::tidy(get_dd_comparisons(fit))
    cmb <- utils::combn(nrow(emm), 2L)
    expected <- (emm$k_log[cmb[1L, ]] - emm$k_log[cmb[2L, ]]) / log(10)
    expect_equal(td$estimate, expected, tolerance = 1e-8)
  })

  it("exponentiate = TRUE returns base-invariant ratios with std.error NA", {
    set.seed(33)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 90, family = "gaussian", equation = "mazur",
      n_conditions = 3, delta_k = 0.6, seed = 33
    )
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = "condition", multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit)
    td <- generics::tidy(res)
    tde <- generics::tidy(res, exponentiate = TRUE)
    expect_equal(tde$estimate, 10^td$estimate, tolerance = 1e-10)
    expect_true(all(is.na(tde$std.error)))
  })

  it("by-columns appear before param when contrast_by is active", {
    set.seed(34)
    dat <- .simulate_dd_ip_mixed(
      n_subjects = 120, family = "gaussian", equation = "mazur",
      n_conditions = 2, delta_k = 0.7, seed = 34
    )
    dat$grp <- factor(rep(c("A", "B"), length.out = nrow(dat)))
    fit <- fit_dd_tmb(dat, equation = "mazur", family = "gaussian",
                      factors = c("condition", "grp"), factor_interaction = TRUE,
                      multi_start = FALSE, verbose = 0)
    res <- get_dd_comparisons(fit, compare_specs = ~ condition + grp,
                              contrast_by = "grp")
    td <- generics::tidy(res)
    expect_true("grp" %in% names(td))
    expect_lt(match("grp", names(td)), match("param", names(td)))
  })
})
```

- [ ] **Step 2: Run the test, expect failure (no tidy method).**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: `no applicable method for 'tidy' applied to an object of class "beezdiscounting_comparison"`.

- [ ] **Step 3: Implement the flattener + `tidy.beezdiscounting_comparison()` in `R/dd-comparisons.R`.**

  Single-parameter port of beezdemand `.beezdemand_comparison_flat()` +
  `tidy.beezdemand_comparison()` (`R/mixed-methods.R:1793-1913`). Adaptations:
  backend is always `"tmb"` (drop the NLME branch and `cb_map` defensiveness);
  the comparison is always single-keyed (`"k"`), but the loop over `names(x)`
  is kept so the schema generalizes; by-column source name is the effective
  factor name directly (no collapse map). Uses `attr(cl, "std_labels")` for the
  `contrast` column.

```r
# R/dd-comparisons.R  (append)

# TRUE when no contrast_by grouping is active (the "NULL" sentinel string).
.dd_contrast_by_inactive <- function(by_used) {
  is.null(by_used) || identical(by_used, "NULL") || !nzchar(by_used)
}

# Flatten a beezdiscounting_comparison to the cross-backend long schema.
# Single-parameter, TMB-only port of .beezdemand_comparison_flat().
.dd_comparison_flat <- function(x, exponentiate = FALSE) {
  by_used <- attr(x, "contrast_by_used")
  by_active <- !.dd_contrast_by_inactive(by_used)
  by_names <- if (by_active) trimws(strsplit(by_used, ",")[[1]]) else character(0)

  base_cols <- list(
    param = character(), contrast = character(), estimate = numeric(),
    std.error = numeric(), statistic = numeric(), df = numeric(),
    conf.low = numeric(), conf.high = numeric(), p.value = numeric()
  )

  rows <- lapply(names(x), function(p) {
    cl <- x[[p]]$contrasts_log10
    if (is.null(cl) || nrow(cl) == 0L || !("estimate" %in% names(cl))) {
      return(NULL)
    }
    lab <- attr(cl, "std_labels")
    if (is.null(lab) || length(lab) != nrow(cl)) lab <- cl$contrast
    base <- tibble::tibble(
      param = p, contrast = lab,
      estimate = cl$estimate, std.error = cl$std.error,
      statistic = cl$statistic, df = cl$df,
      conf.low = cl$conf.low, conf.high = cl$conf.high,
      p.value = cl$p.value
    )
    if (by_active) {
      by_cols <- lapply(by_names, function(nm) {
        if (nm %in% names(cl)) as.character(cl[[nm]]) else rep(NA_character_, nrow(cl))
      })
      base <- dplyr::bind_cols(
        tibble::as_tibble(stats::setNames(by_cols, by_names)), base
      )
    }
    base
  })

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0L) {
    out <- tibble::as_tibble(base_cols)
    if (by_active) {
      by_empty <- stats::setNames(rep(list(character()), length(by_names)), by_names)
      out <- dplyr::bind_cols(tibble::as_tibble(by_empty), out)
    }
  }

  if (isTRUE(exponentiate)) {
    out$estimate <- 10^out$estimate
    out$conf.low <- 10^out$conf.low
    out$conf.high <- 10^out$conf.high
    out$std.error <- NA_real_
  }
  out
}

#' Tidy a discounting comparison into a flat contrasts frame
#'
#' @description
#' [generics::tidy()] method for \code{beezdiscounting_comparison} objects
#' returned by \code{\link{get_dd_comparisons}}. Returns a flat long tibble:
#' one row per contrast, on the \code{log10 k} scale (or base-invariant ratios
#' when \code{exponentiate = TRUE}).
#'
#' @param x A \code{beezdiscounting_comparison} object.
#' @param exponentiate Logical. If \code{TRUE}, return ratios
#'   (\code{estimate = 10^estimate}, CIs back-transformed); \code{std.error}
#'   becomes \code{NA} (broom's exponentiated-fit convention). Default
#'   \code{FALSE}.
#' @param ... Unused.
#'
#' @return A tibble with columns \code{param} (always \code{"k"}),
#'   \code{contrast}, \code{estimate}, \code{std.error}, \code{statistic}
#'   (asymptotic \emph{z}), \code{df} (\code{Inf}), \code{conf.low},
#'   \code{conf.high}, \code{p.value}. When \code{contrast_by} was active,
#'   leading by-column(s) precede \code{param}.
#'
#' @importFrom generics tidy
#' @export
tidy.beezdiscounting_comparison <- function(x, exponentiate = FALSE, ...) {
  .dd_comparison_flat(x, exponentiate = exponentiate)
}
```

  Ensure `generics` is imported. Add to `DESCRIPTION` `Imports:` if absent:
  `generics`. The `#' @importFrom generics tidy` above re-exports the generic;
  also add a standalone re-export so `tidy()` is callable without
  `generics::`:

```r
# R/dd-comparisons.R  (append)
#' @importFrom generics tidy
#' @export
generics::tidy
```

- [ ] **Step 4: Run the test, expect pass.**

  Run: `Rscript -e 'devtools::load_all("/Users/brent/Dropbox/GIT/beezdiscounting"); testthat::test_file("/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-comparisons.R")'`

  Expected: all `tidy.beezdiscounting_comparison()` `it()` blocks pass, and the full `test-dd-comparisons.R` file is green.

- [ ] **Step 5: Regenerate docs/NAMESPACE, lint, and commit.**

  Run:
  ```
  Rscript -e 'devtools::document("/Users/brent/Dropbox/GIT/beezdiscounting")'
  Rscript -e 'lintr::lint("/Users/brent/Dropbox/GIT/beezdiscounting/R/dd-comparisons.R")'
  git -C /Users/brent/Dropbox/GIT/beezdiscounting add R/dd-comparisons.R tests/testthat/test-dd-comparisons.R NAMESPACE DESCRIPTION man/tidy.beezdiscounting_comparison.Rd man/reexports.Rd
  git -C /Users/brent/Dropbox/GIT/beezdiscounting commit -m "feat(emmeans): add tidy.beezdiscounting_comparison flat frame

Backend-agnostic flattener (param='k') with structured std_labels,
contrast_by columns before param, and exponentiate ratios; re-export the
generics::tidy generic.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
  ```
## Phase S — Simulation, tie-out, vignette, CI, release

> Surface-level deliverables that sit on top of the already-built TMB engine
> (`fit_dd_tmb`, the S3 methods, `get_dd_param_emms`/`get_dd_comparisons`) and the validator
> (`.dd_validate_ip`). This phase adds the data-generating simulator, the NLS/two-stage
> tie-out tests, backfilled coverage for the legacy NLS path, the shipped vignette + bundled
> Jarvis fixture, the TMB-aware CI workflow, and the release bookkeeping (NEWS, version 0.4.0,
> NAMESPACE regen).
>
> All compiled-template tests open with `skip_on_cran()` and `skip_if_not_installed("TMB")`.
> Tests use testthat 3e BDD (`describe()`/`it()`). Every commit uses targeted `git add`
> (never `-A`) and ends the message with the pinned `Co-Authored-By` trailer.

### Task S.1: `.simulate_dd_ip_mixed()` — SLT inverse-CDF + Gaussian DGP

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/simulate-dd-mixed.R`
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-simulate-dd-mixed.R`
- Modify (helper, append): `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/helper-dd-sim.R`

- [ ] **Step 1: Write the failing test** for the column contract, factor branch, and the SLT moment of the inverse-CDF draw. Create `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-simulate-dd-mixed.R`:

```r
describe(".simulate_dd_ip_mixed()", {
  it("returns the long-format column contract (id, x, y) with no condition by default", {
    sim <- .simulate_dd_ip_mixed(n_subjects = 20, seed = 1)
    expect_s3_class(sim, "tbl_df")
    expect_identical(names(sim), c("id", "x", "y"))
    expect_s3_class(sim$id, "factor")
    expect_equal(nlevels(sim$id), 20L)
    expect_true(all(sim$y >= 0 & sim$y <= 1))
    # 7 default delays per subject
    expect_equal(nrow(sim), 20L * 7L)
  })

  it("adds a condition factor when n_conditions > 1", {
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 12, n_conditions = 2, delta_k = c(0, log(3)), seed = 2
    )
    expect_identical(names(sim), c("id", "condition", "x", "y"))
    expect_s3_class(sim$condition, "factor")
    expect_equal(nlevels(sim$condition), 2L)
  })

  it("errors when delta_k length does not match n_conditions", {
    expect_error(
      .simulate_dd_ip_mixed(n_subjects = 5, n_conditions = 3, delta_k = c(0, 1)),
      "delta_k"
    )
  })

  it("SLT draws recover the population mean curve in expectation (mazur)", {
    set.seed(99)
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 400, log_k_pop = log(0.01), sigma_u = 1e-6, phi = 40,
      family = "sltb", equation = "mazur", seed = 99
    )
    # with sigma_u ~ 0 every subject shares k = exp(log_k_pop); mean y at each delay
    # should track mu = 1/(1 + k*x) to within Monte-Carlo error
    k <- 0.01
    by_delay <- tapply(sim$y, sim$x, mean)
    mu_true <- 1 / (1 + k * as.numeric(names(by_delay)))
    expect_lt(max(abs(by_delay - mu_true)), 0.02)
  })

  it("gaussian draws clamp to [0,1]", {
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 50, family = "gaussian", sigma_e = 0.4, seed = 7
    )
    expect_true(all(sim$y >= 0 & sim$y <= 1))
  })
})
```

- [ ] **Step 2: Run the test (expect fail — function does not exist).**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-simulate-dd-mixed.R")'
```
Expected output: errors `could not find function ".simulate_dd_ip_mixed"` (all `it()` blocks fail).

- [ ] **Step 3: Implement `.simulate_dd_ip_mixed()`.** Create `/Users/brent/Dropbox/GIT/beezdiscounting/R/simulate-dd-mixed.R` (complete; SLT inverse-CDF exactly as in `dev/sltb-verification/verify_sltb.R::rslt` and the contract):

```r
#' Simulate IP-family mixed-effects discounting data
#'
#' Generates long-format indifference-point data from the mixed-effects
#' discounting model used by [fit_dd_tmb()]. Each subject `i` has a random
#' discount rate `log k_i = log_k_pop + delta_k[condition] + u_i`,
#' `u_i ~ N(0, sigma_u^2)`. The mean indifference proportion at delay `x` is the
#' discounting function `mu` (Mazur hyperbola or exponential), and observed `y`
#' is drawn from the scale-location-truncated beta (`family = "sltb"`) via the
#' inverse-CDF on the truncated beta, or from a clamped Gaussian
#' (`family = "gaussian"`).
#'
#' The SLT draw uses the same constants as the C++ template and the verified
#' reference density: `s = 1.0000001`, `l = 1e-8`, with
#' `a = mu * phi`, `b = (1 - mu) * phi`, and
#' `y = (qbeta(U, a, b) - l) * s` where
#' `U ~ Uniform(pbeta(l, a, b), pbeta(1/s + l, a, b))`.
#'
#' @param n_subjects Integer; number of subjects.
#' @param delays Numeric vector of delays (days) each subject is observed at.
#' @param log_k_pop Numeric; population intercept on `log k`.
#' @param sigma_u Numeric; SD of the subject random intercept on `log k`.
#' @param phi Numeric; SLT-beta precision (`family = "sltb"`).
#' @param sigma_e Numeric; residual SD on `y` (`family = "gaussian"`).
#' @param family One of `"sltb"` (default) or `"gaussian"`.
#' @param equation One of `"mazur"` (default) or `"exponential"`.
#' @param n_conditions Integer; number of between-subject condition levels. When
#'   `> 1`, a `condition` factor is added and subjects are split across levels.
#' @param delta_k Numeric vector of length `n_conditions`; per-condition shift on
#'   `log k` (the first element is typically `0` for the reference level).
#'   Required (non-`NULL`) when `n_conditions > 1`.
#' @param seed Optional integer seed.
#'
#' @return A [tibble][tibble::tibble] with columns `id` (factor),
#'   `condition` (factor; only when `n_conditions > 1`), `x` (delay), and
#'   `y` (indifference proportion in `[0, 1]`).
#'
#' @keywords internal
#' @importFrom stats rnorm runif qbeta pbeta
#' @importFrom tibble tibble
.simulate_dd_ip_mixed <- function(
  n_subjects = 60,
  delays = c(7, 30, 180, 365, 730, 1460, 2920),
  log_k_pop = log(0.01),
  sigma_u = 0.6,
  phi = 10,
  sigma_e = 0.1,
  family = c("sltb", "gaussian"),
  equation = c("mazur", "exponential"),
  n_conditions = 1,
  delta_k = NULL,
  seed = NULL
) {
  family <- match.arg(family)
  equation <- match.arg(equation)
  if (!is.null(seed)) set.seed(seed)

  if (n_conditions > 1 && is.null(delta_k)) {
    stop("`delta_k` must be supplied (length `n_conditions`) when `n_conditions > 1`.",
         call. = FALSE)
  }
  if (is.null(delta_k)) delta_k <- rep(0, n_conditions)
  if (length(delta_k) != n_conditions) {
    stop("`delta_k` must have length `n_conditions`.", call. = FALSE)
  }

  s <- 1.0000001
  l <- 1e-8

  # assign each subject to a condition (balanced round-robin)
  cond_idx <- rep_len(seq_len(n_conditions), n_subjects)
  cond_lab <- factor(paste0("C", cond_idx), levels = paste0("C", seq_len(n_conditions)))

  u <- stats::rnorm(n_subjects, 0, sigma_u)
  log_k <- log_k_pop + delta_k[cond_idx] + u
  k <- exp(log_k)

  n_delays <- length(delays)
  id <- factor(rep(seq_len(n_subjects), each = n_delays))
  x <- rep(delays, times = n_subjects)
  k_long <- rep(k, each = n_delays)

  mu <- if (equation == "mazur") {
    1 / (1 + k_long * x)
  } else {
    exp(-k_long * x)
  }
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)

  if (family == "sltb") {
    a <- mu * phi
    b <- (1 - mu) * phi
    lo <- stats::pbeta(l, a, b)
    hi <- stats::pbeta(1 / s + l, a, b)
    uu <- stats::runif(length(mu), lo, hi)
    y <- (stats::qbeta(uu, a, b) - l) * s
  } else {
    y <- stats::rnorm(length(mu), mu, sigma_e)
  }
  y <- pmin(pmax(y, 0), 1)

  if (n_conditions > 1) {
    tibble::tibble(
      id = id,
      condition = rep(cond_lab, each = n_delays),
      x = x,
      y = y
    )
  } else {
    tibble::tibble(id = id, x = x, y = y)
  }
}
```

- [ ] **Step 4: Add a reusable simulate helper** so the recovery/tie-out tests share one constructor. Append to `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/helper-dd-sim.R` (create the file if it does not exist yet from earlier phases):

```r
# Deterministic small SLT-beta mixed sim used across recovery + tie-out tests.
.dd_sim_fixture <- function(family = "sltb", equation = "mazur",
                            n_subjects = 80, seed = 20260606) {
  .simulate_dd_ip_mixed(
    n_subjects = n_subjects,
    log_k_pop = log(0.01),
    sigma_u = 0.6,
    phi = 12,
    sigma_e = 0.1,
    family = family,
    equation = equation,
    seed = seed
  )
}
```

- [ ] **Step 5: Add the known-truth recovery test through `fit_dd_tmb()`.** Append to `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-simulate-dd-mixed.R`:

```r
describe(".simulate_dd_ip_mixed() recovery through fit_dd_tmb()", {
  it("recovers population k within 0.15 relative (sltb, mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- .dd_sim_fixture(family = "sltb", equation = "mazur", seed = 101)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    k_hat <- exp(unname(fit$model$coefficients[["beta_k"]]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.15)
  })

  it("recovers population k within 0.15 relative (gaussian, mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- .dd_sim_fixture(family = "gaussian", equation = "mazur", seed = 102)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "gaussian",
                      multi_start = TRUE, verbose = 0)
    k_hat <- exp(unname(fit$model$coefficients[["beta_k"]]))
    expect_lt(abs(k_hat - 0.01) / 0.01, 0.15)
  })

  it("recovers an exponential rate by rank correlation over a k grid", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    # exponential is less identifiable; use the identifiable regime + looser bar
    sim <- .simulate_dd_ip_mixed(
      n_subjects = 120, log_k_pop = log(3e-4), sigma_u = 0.6, phi = 12,
      family = "sltb", equation = "exponential", seed = 103
    )
    fit <- fit_dd_tmb(sim, equation = "exponential", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    sp <- fit$subject_pars
    # per-subject two-stage k for the same subjects, ranked agreement
    ts <- vapply(split(sim, sim$id), function(d) {
      m <- tryCatch(
        stats::nls(y ~ exp(-k * x), data = d, start = list(k = 3e-4)),
        error = function(e) NULL
      )
      if (is.null(m)) NA_real_ else unname(coef(m)[["k"]])
    }, numeric(1))
    ok <- is.finite(ts) & ts > 0
    expect_gt(cor(log(sp$k[ok]), log(ts[ok])), 0.93)
  })
})
```

Note: `beta_k` is the FE name from the contract; for the intercept-only design `model.matrix` yields a single `(Intercept)` column whose coefficient is named `beta_k` in `fit$model$coefficients` per the TMB phase. If that phase names the single intercept `"(Intercept)"` instead, change `coefficients[["beta_k"]]` to `coefficients[["(Intercept)"]]` — confirm against the symbol table the TMB phase exports and keep one form.

- [ ] **Step 6: Run the full file (expect pass).**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-simulate-dd-mixed.R")'
```
Expected output: all `it()` blocks pass (compiled-template recovery blocks pass once `TMB` is installed; otherwise they skip, never fail).

- [ ] **Step 7: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add R/simulate-dd-mixed.R tests/testthat/test-simulate-dd-mixed.R tests/testthat/helper-dd-sim.R && \
git commit -m "feat(sim): add .simulate_dd_ip_mixed SLT/Gaussian DGP with recovery tests

SLT draws via inverse-CDF on the truncated beta (s=1.0000001, l=1e-8);
Gaussian draws clamped to [0,1]; optional between-subject condition factor via
delta_k. Recovery of population k within 0.15 relative through fit_dd_tmb for
both families (mazur) and rank correlation > 0.93 for the exponential.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task S.2: `test-dd-nls-tieout.R` — Gaussian ties to per-subject NLS; SLTB intercept ≈ median(two-stage k)

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-nls-tieout.R`

- [ ] **Step 1: Write the tie-out test.** The key statistical point: the population FE recovers a *central* k on the log scale, so `exp(beta_k_intercept)` is a **geometric mean / median** of subject k, NOT the arithmetic mean — under log-normal RE the two differ by `exp(sigma^2/2)`. Create `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-dd-nls-tieout.R`:

```r
# Per-subject two-stage k via NLS (multi-start; single start is fragile).
.two_stage_nls_k <- function(data, equation = "mazur") {
  fo <- if (equation == "mazur") y ~ 1 / (1 + k * x) else y ~ exp(-k * x)
  vapply(split(data, data$id), function(d) {
    for (st in c(exp(-10), 1e-4, 1e-3, 1e-2, 1e-1)) {
      m <- tryCatch(
        stats::nls(fo, data = d, start = list(k = st)),
        error = function(e) NULL
      )
      if (!is.null(m)) return(unname(coef(m)[["k"]]))
    }
    NA_real_
  }, numeric(1))
}

describe("fit_dd_tmb() ties out to per-subject NLS", {
  it("gaussian family agrees with two-stage NLS k on the log scale (mazur)", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- .dd_sim_fixture(family = "gaussian", equation = "mazur", seed = 2001)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "gaussian",
                      multi_start = TRUE, verbose = 0)
    sp <- fit$subject_pars                       # per-subject k (BLUP-shrunken)
    ts <- .two_stage_nls_k(sim, equation = "mazur")
    ok <- is.finite(ts) & ts > 0 & is.finite(sp$k) & sp$k > 0
    # gaussian + same mean function => subject k matches NLS k closely in rank/level
    expect_gt(cor(log(sp$k[ok]), log(ts[ok])), 0.95)
  })

  it("population intercept is the GEOMETRIC mean (median), not the arithmetic mean", {
    skip_on_cran()
    skip_if_not_installed("TMB")
    sim <- .dd_sim_fixture(family = "sltb", equation = "mazur", seed = 2002)
    fit <- fit_dd_tmb(sim, equation = "mazur", family = "sltb",
                      multi_start = TRUE, verbose = 0)
    k_pop <- exp(unname(fit$model$coefficients[["beta_k"]]))
    ts <- .two_stage_nls_k(sim, equation = "mazur")
    ts <- ts[is.finite(ts) & ts > 0]

    geo  <- exp(mean(log(ts)))     # geometric mean of two-stage k
    med  <- median(ts)            # median (= geometric mean under symmetric log RE)
    arith <- mean(ts)             # arithmetic mean (inflated by exp(sigma^2/2))

    # exp(beta0) tracks the geometric mean / median, within 0.20 relative
    expect_lt(abs(k_pop - geo) / geo, 0.20)
    expect_lt(abs(k_pop - med) / med, 0.25)

    # and it is demonstrably BELOW the arithmetic mean by ~ exp(sigma^2/2):
    # with sigma_u = 0.6, exp(sigma^2/2) = exp(0.18) ~ 1.197, so arith > geo.
    expect_gt(arith, geo)
    sigma_u_hat <- unname(fit$model$coefficients[["sigma_u"]])
    expect_lt(abs((arith / geo) - exp(sigma_u_hat^2 / 2)) / exp(sigma_u_hat^2 / 2), 0.35)
  })
})
```

Note: `sigma_u` appears in `fit$model$coefficients` as the ADREPORTed `sigma_u = exp(log_sigma_u)` per the TMB phase. If that phase exposes it only under `fit$model$coefficients[["log_sigma_u"]]`, replace `sigma_u_hat <- ...["sigma_u"]` with `sigma_u_hat <- exp(...[["log_sigma_u"]])`. Reconcile against the TMB phase symbol table; keep exactly one form.

- [ ] **Step 2: Run (expect pass once the TMB engine exists; skip without TMB).**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-dd-nls-tieout.R")'
```
Expected output: both `it()` blocks pass (or skip cleanly if `TMB` not installed). The geometric-mean block confirms `exp(beta0)` < arithmetic mean and matches `exp(sigma_u^2/2)`.

- [ ] **Step 3: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add tests/testthat/test-dd-nls-tieout.R && \
git commit -m "test(tieout): tie fit_dd_tmb to per-subject NLS via geometric-mean k

Gaussian family agrees with two-stage NLS k on the log scale; the SLTB
population intercept exp(beta0) is the geometric mean / median of the two-stage
k (NOT the arithmetic mean, which is inflated by exp(sigma^2/2) under the
log-normal random intercept).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task S.3: Backfill `test-fit_dd.R` for the legacy NLS path (`fit_dd`/`results_dd`)

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd.R`
- Reference (read-only, current behavior): `/Users/brent/Dropbox/GIT/beezdiscounting/R/fitting.R:29` (`fit_dd`), `:248` (`results_dd`)

Current behavior to lock in (from `R/fitting.R`):
- `fit_dd(dat, equation, method)` returns a length-3 list of class `c("fit_dd", "list")`: `[[1]]` the fit(s), `[[2]]` the data (with `id` coerced to factor), `[[3]]` the method string.
- `method = "pooled"`/`"agg"` and `"mean"` fit a single `minpack.lm::nlsLM` (`[[1]]` is one model); `"ts"`/`"two stage"` returns a named list (one `purrr::safely` result per `id`, each `$result`/`$error`).
- `results_dd()` for pooled/mean returns a tibble with `method`, `term`, `estimate`, `std.error`, `statistic`, `p.value`, `R2`, `conf_low`, `conf_high` (mean additionally adds `auc*` columns); for two-stage returns one row per `id` with `method`, `id`, `term`, ... , `R2`, `auc*`, `conf_low`, `conf_high`.
- Equation accepts mixed case (`"mazur"/"Mazur"/"hyperbolic"/"Hyperbolic"/"exponential"/"Exponential"`); method accepts mixed case too. Out-of-domain `equation`/`method` trip `stopifnot`.

- [ ] **Step 1: Write the characterization tests.** Create `/Users/brent/Dropbox/GIT/beezdiscounting/tests/testthat/test-fit_dd.R`:

```r
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
```

- [ ] **Step 2: Run (expect pass — these characterize existing behavior).**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-fit_dd.R")'
```
Expected output: all `it()` blocks pass against the current `fitting.R`. If any column name assertion fails, fix the *test* to match the actual `results_dd()` output (this is a characterization backfill of existing behavior, not a behavior change) — read `R/fitting.R:248` and adjust the asserted column set.

- [ ] **Step 3: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add tests/testthat/test-fit_dd.R && \
git commit -m "test(fit_dd): backfill characterization tests for the legacy NLS path

Locks in fit_dd() (pooled/mean/two-stage shapes, mixed-case args, stopifnot
guards) and results_dd() (tidy columns, AUC for mean, one row per id for
two-stage) before the TMB tier ships alongside it.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task S.4: Bundled Jarvis fixture (`data-raw/jarvis2019.R` → `data/jarvis2019.rda`) + roxygen doc

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/data-raw/jarvis2019.R`
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/R/jarvis2019.R` (roxygen dataset doc)
- Create (build artifact): `/Users/brent/Dropbox/GIT/beezdiscounting/data/jarvis2019.rda`

The source CSV (`Data_brent.csv`) is **not in the repo** (lives at `~/Downloads/Data_brent.csv`). `PID` is non-unique (78 unique / 126 rows), so the fixture **must use a unique synthetic subject id** (row index) and document that each row is one subject-record. The build script is documented and idempotent; if the CSV is absent it errors with a clear message rather than producing a silent empty fixture.

- [ ] **Step 1: Write the data-raw build script.** Create `/Users/brent/Dropbox/GIT/beezdiscounting/data-raw/jarvis2019.R`:

```r
## ----------------------------------------------------------------------------
## Build the bundled, de-duplicated Jarvis (2019) delay-discounting fixture.
##
## SOURCE: Data_brent.csv (Jarvis et al. 2019; N = 126, delays 1wk..8yr).
## The CSV is NOT shipped in the package (private source data); it lives at
## ~/Downloads/Data_brent.csv. Run this script manually to (re)build
## data/jarvis2019.rda. PID is NON-UNIQUE (78 unique / 126 rows), so we assign
## a synthetic unique subject id (the row index) and treat each row as one
## subject-record. A mixed model that grouped duplicated PIDs as one subject
## would mis-group; the unique `id` here prevents that.
## ----------------------------------------------------------------------------

csv_path <- path.expand("~/Downloads/Data_brent.csv")
if (!file.exists(csv_path)) {
  stop(
    "Source CSV not found at ", csv_path, ".\n",
    "Place the Jarvis (2019) Data_brent.csv there and re-run this script.\n",
    "The CSV is not distributed with the package.",
    call. = FALSE
  )
}

d <- utils::read.csv(csv_path, check.names = FALSE)
delays <- c(7, 30, 180, 365, 730, 1460, 2920)  # 1wk,1mo,6mo,1yr,2yr,4yr,8yr (days)

# Each row is one subject-record; PID is non-unique so we use the row index.
n <- nrow(d)
ip_cols <- d[, -1, drop = FALSE]               # drop the PID column
stopifnot(ncol(ip_cols) == length(delays))

jarvis2019 <- tibble::tibble(
  id = factor(rep(seq_len(n), each = length(delays))),
  x  = rep(delays, times = n),
  y  = as.numeric(t(as.matrix(ip_cols)))
)

# Sanity: boundary IPs are the point of this fixture.
message(sprintf(
  "jarvis2019: %d subjects, %d obs, IP in [%.3f, %.3f]; #(y==0)=%d, #(y==1)=%d",
  n, nrow(jarvis2019), min(jarvis2019$y), max(jarvis2019$y),
  sum(jarvis2019$y == 0), sum(jarvis2019$y == 1)
))

usethis::use_data(jarvis2019, overwrite = TRUE)
```

- [ ] **Step 2: Write the dataset roxygen doc.** Create `/Users/brent/Dropbox/GIT/beezdiscounting/R/jarvis2019.R`:

```r
#' Jarvis (2019) delay-discounting indifference points
#'
#' Long-format delay-discounting indifference-point data from Jarvis et al.
#' (2019), used to demonstrate the scale-location-truncated beta (SLT-beta)
#' family's handling of boundary indifference points (exactly 0 or 1). Many
#' observations sit at the bounds, where ordinary beta regression is undefined.
#'
#' The original source has a non-unique participant identifier (78 unique values
#' across 126 records); here each record is treated as one subject-record and is
#' assigned a unique synthetic `id` (the row index of the source file). Use this
#' `id` as the subject grouping for [fit_dd_tmb()].
#'
#' @format A [tibble][tibble::tibble] with 882 rows (126 subjects x 7 delays)
#'   and 3 columns:
#' \describe{
#'   \item{id}{Factor; unique subject-record identifier.}
#'   \item{x}{Numeric; delay in days (7, 30, 180, 365, 730, 1460, 2920 =
#'     1 week through 8 years).}
#'   \item{y}{Numeric; indifference point as a proportion of the larger-later
#'     reward, in the closed interval `[0, 1]` (0 and 1 occur and are valid).}
#' }
#'
#' @source Jarvis, B. P., et al. (2019). Delay-discounting dataset
#'   (N = 126; 1 week - 8 years). Source CSV not redistributed with the package;
#'   see `data-raw/jarvis2019.R`.
#' @keywords datasets
"jarvis2019"
```

- [ ] **Step 3: Build the fixture (manual; requires the CSV).**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript data-raw/jarvis2019.R
```
Expected output: a `jarvis2019: 126 subjects, 882 obs, IP in [0.000, 1.000]; #(y==0)=..., #(y==1)=...` message and `data/jarvis2019.rda` written. If the CSV is absent, the script stops with the documented "Source CSV not found" error — in that case the fixture/vignette demo cannot be rebuilt and the vignette boundary chunks must remain guarded (see Task S.5 below). Record this state for the integrator.

- [ ] **Step 4: Verify the fixture loads and has the boundary property.**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'load("data/jarvis2019.rda"); stopifnot(identical(names(jarvis2019), c("id","x","y"))); cat("subjects:", nlevels(jarvis2019$id), "y0:", sum(jarvis2019$y==0), "y1:", sum(jarvis2019$y==1), "\n")'
```
Expected output: `subjects: 126 y0: <n> y1: <n>` with at least one of `y0`/`y1` positive (boundary IPs present).

- [ ] **Step 5: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add data-raw/jarvis2019.R R/jarvis2019.R data/jarvis2019.rda && \
git commit -m "data(jarvis2019): add de-duplicated Jarvis (2019) boundary fixture

data-raw/jarvis2019.R builds data/jarvis2019.rda from the (non-shipped)
Data_brent.csv, assigning a unique synthetic subject id (row index) because the
source PID is non-unique. Long format id/x/y in [0,1] with boundary IPs (0 and
1) for the SLT-beta demonstration. Roxygen dataset doc added.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task S.5: Convert the draft into `vignettes/sltb-discounting.Rmd` (flip API chunks to `eval=TRUE`)

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/vignettes/sltb-discounting.Rmd`
- Reference (read-only, source): `/Users/brent/Dropbox/GIT/beezdiscounting/dev/sltb-verification/vignette-draft.Rmd`

Adaptations from the draft:
1. **Data loading** — replace the `~/Downloads/Data_brent.csv` `if (file.exists(...))` block with `data(jarvis2019)`; rename `long`→`jarvis2019` and `delay`→`x`, `IP`→`y` throughout the boundary/fit chunks to match the bundled fixture columns.
2. **API chunks** (`fit-tmb`, `emmeans`) — flip `eval = FALSE` → `eval = TRUE` (the engine now exists). Use `equation`/`family`/`random_effects`/`factors` exactly per the pinned signature.
3. Keep the `slt_logpdf`/`beta_logpdf` helper definitions inline (the vignette is self-contained and must not reach into non-exported internals).
4. Guard the **per-subject** comparison chunk and the API chunks behind `requireNamespace("TMB", quietly = TRUE)` so the vignette still builds on machines without a compiler/TMB (CRAN safe).

- [ ] **Step 1: Write the vignette.** Create `/Users/brent/Dropbox/GIT/beezdiscounting/vignettes/sltb-discounting.Rmd` (full content; sections 1–3, 6 run on the bundled fixture, 4–5 run when TMB is available):

````markdown
---
title: "Modeling delay-discounting indifference points with bounded error distributions"
author: "Brent Kaplan"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Modeling delay-discounting indifference points with bounded error distributions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.2)
set.seed(1)
has_tmb <- requireNamespace("TMB", quietly = TRUE)
library(beezdiscounting)
```

## Why indifference points need a bounded error distribution

A delay-discounting task asks a person to trade a smaller-sooner reward against a
larger-later one across a range of delays. The **indifference point** (IP) at each delay is
the smaller-sooner amount that feels equivalent to the later reward, expressed as a
*proportion of that later reward*. By construction an IP lives on the closed interval
$[0, 1]$:

- $y = 1$ — the person is indifferent only when the immediate amount equals the full later
  reward: **no discounting**.
- $y = 0$ — any immediate amount, however small, is preferred: **complete discounting**.

Real data pile up at both ends. In the Jarvis et al. (2019) data shipped with this package
(`data(jarvis2019)`), a meaningful fraction of observations sit *exactly* at 0 or 1. That
single fact rules out the two error distributions people reach for first:

- **Gaussian / nonlinear least squares** treats the residuals as unbounded and
  homoscedastic. It will happily predict an IP below 0 or above 1, and it assumes the noise
  is the same size near the bounds as in the middle — which is never true for a proportion.
- **Standard beta regression** is the natural model for a proportion, but its density is
  *undefined* at exactly 0 and 1 (the log-likelihood is $-\infty$). Every boundary
  observation has to be discarded or nudged, which biases the fit.

The **Scale-Location-Truncated beta (SLT-beta)** distribution (Kim, Koffarnus & Franck, 2024;
Kim, Kaplan, Koffarnus & Franck, 2025) keeps the beta's natural fit to proportion data but
assigns *finite probability* to the endpoints, so 0 and 1 are modeled rather than thrown away.

## What the SLT-beta density is

Write the mean of the IP at delay $D$ as a discounting function of a single rate $k$ — Mazur's
hyperbola $\mu = 1/(1 + kD)$ or the exponential $\mu = e^{-kD}$ — and let $\phi$ be a
precision (larger $\phi$ = tighter around the curve). With shape parameters
$a = \mu\phi$ and $b = (1-\mu)\phi$, and two *fixed* micro-constants $s = 1{+}10^{-7}$ and
$l = 10^{-8}$, the log-density of an observed IP $y \in [0,1]$ is

$$
\log f(y \mid \mu, \phi) =
\underbrace{\log\Gamma(a{+}b) - \log\Gamma(a) - \log\Gamma(b)}_{-\log B(a,b)}
+ (a{-}1)\log\!\big(\tfrac{y}{s}{+}l\big)
+ (b{-}1)\log\!\big(1 - \tfrac{y}{s}{+}l\big)
- \log s
- \log Z,
$$

where the truncation normalizer is $Z = \mathrm{pbeta}(\tfrac{1}{s}{+}l,\,a,\,b) -
\mathrm{pbeta}(l,\,a,\,b)$.

The intuition is a *microscopic stretch*: the data window $y\in[0,1]$ is mapped just inside
the open interval where the ordinary beta density is finite, so the kernel is never evaluated
at exactly 0 or 1. As $s\to 1$ and $l\to 0$ the SLT-beta collapses to the ordinary beta.

```{r density, fig.cap = "SLT-beta (solid) vs. ordinary beta (dashed). The SLT-beta carries finite density to y = 0 and y = 1; the ordinary beta diverges or vanishes there."}
slt_logpdf <- function(y, mu, phi, s = 1.0000001, l = 1e-8) {
  a <- mu * phi; b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(a) - lgamma(b) +
    (a - 1) * log(y / s + l) + (b - 1) * log(1 - (y / s + l)) -
    log(s) - log(pbeta(1 / s + l, a, b) - pbeta(l, a, b))
}
yy <- seq(0, 1, length.out = 400)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (p in list(c(mu = 0.5, phi = 4), c(mu = 0.15, phi = 6))) {
  plot(yy, exp(slt_logpdf(yy, p["mu"], p["phi"])), type = "l", lwd = 2,
       xlab = "indifference point", ylab = "density",
       main = sprintf("mu = %.2f, phi = %g", p["mu"], p["phi"]))
  lines(yy, dbeta(yy, p["mu"] * p["phi"], (1 - p["mu"]) * p["phi"]), lty = 2)
  points(c(0, 1), exp(slt_logpdf(c(0, 1), p["mu"], p["phi"])), pch = 16)
}
par(op)
```

Because the variance of a beta is $\mu(1-\mu)/(\phi+1)$, the SLT-beta automatically gives
**less** noise near the floor and ceiling and **more** in the middle — the heteroscedastic
pattern real IP data actually show, and exactly what a Gaussian model gets wrong.

## The boundary problem, on real data

```{r load-data}
data(jarvis2019)
knitr::kable(as.data.frame(table(IP = jarvis2019$y)),
             caption = "Observed indifference points")
```

The titration grid produces many observations at exactly 0 and 1. Counting subjects with at
least one boundary value, and asking which error models can even evaluate them:

```{r boundary-demo}
beta_logpdf <- function(y, mu, phi) {          # ordinary beta: -Inf at y in {0,1}
  a <- mu * phi; b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(a) - lgamma(b) + (a - 1) * log(y) + (b - 1) * log(1 - y)
}
ids <- levels(jarvis2019$id)
yof <- function(i) jarvis2019$y[jarvis2019$id == i]
bnd     <- vapply(ids, function(i) any(yof(i) %in% c(0, 1)), logical(1))
beta_ok <- vapply(ids, function(i) all(is.finite(beta_logpdf(yof(i), 0.5, 8))), logical(1))
slt_ok  <- vapply(ids, function(i) all(is.finite(slt_logpdf(yof(i), 0.5, 8))), logical(1))
data.frame(
  subjects                       = length(ids),
  with_boundary_IP               = sum(bnd),
  evaluable_under_ordinary_beta  = sum(beta_ok),
  evaluable_under_SLT_beta       = sum(slt_ok)
)
```

Every subject with a boundary value is lost to ordinary beta regression and retained by the
SLT-beta. Fitting each subject (estimating $\log k$ and $\log\phi$ so both stay positive) and
comparing the recovered rate to nonlinear least squares shows the SLT-beta agrees with NLS on
the bulk of subjects *and* succeeds on the boundary subjects NLS-plus-beta could not model:

```{r per-subject-fit}
fit_slt <- function(y, D) {                      # per-subject SLT MLE (log k, log phi)
  nll <- function(th) {
    mu <- pmin(pmax(1 / (1 + exp(th[1]) * D), 1e-6), 1 - 1e-6)
    v <- -sum(slt_logpdf(y, mu, exp(th[2]))); if (is.finite(v)) v else 1e10
  }
  o <- optim(c(log(0.01), log(8)), nll, control = list(reltol = 1e-11))
  exp(o$par[1])
}
nls_k <- function(y, D) {                        # multi-start NLS (single start is fragile)
  for (st in c(exp(-10), 1e-3, 1e-2, 1e-1)) {
    m <- tryCatch(nls(y ~ 1 / (1 + k * D), start = list(k = st)), error = function(e) NULL)
    if (!is.null(m)) return(coef(m)[["k"]])
  }
  NA_real_
}
xof <- function(i) jarvis2019$x[jarvis2019$id == i]
k_slt <- vapply(ids, function(i) fit_slt(yof(i), xof(i)), numeric(1))
k_nls <- vapply(ids, function(i) nls_k(yof(i), xof(i)), numeric(1))
sprintf("cor(log k): SLT vs NLS = %.3f   |   median k: SLT = %.4f, NLS = %.4f",
        cor(log(k_slt), log(k_nls), use = "complete.obs"),
        median(k_slt), median(k_nls, na.rm = TRUE))
```

A caveat the next section resolves: fit *one subject at a time* and the SLT-beta likelihood
can occasionally run off to a degenerate solution — for a person whose indifference points are
almost all 0s and 1s, a fit with near-zero precision that places all its mass on the endpoints
can score a higher likelihood than any sensible discounting curve. The fix is not to fit
people in isolation.

## Fitting the mixed-effects model

A mixed-effects model shares information across people. The subject random intercept on
$\log k$ acts as a prior that pulls each person's estimate toward the population, so the
degenerate per-subject solutions above are regularized away — one of the main reasons to
prefer `fit_dd_tmb()` over fitting each subject separately.

```{r fit-tmb, eval = has_tmb}
fit <- fit_dd_tmb(
  data           = jarvis2019,       # columns: id, x (delay), y (indifference point)
  equation       = "mazur",          # or "exponential"
  family         = "sltb",           # default; "gaussian" reproduces the NLS-style fit
  random_effects = k ~ 1             # random intercept on log k (fit internally on log scale)
)
summary(fit)
tidy(fit)                            # population log k, precision; back-transformed k
ranef(fit)                           # per-subject discount rates
```

Group contrasts use the same estimated-marginal-means surface as `beezdemand`, computed on
the `log k` scale and back-transformed, so a between-group comparison reads as a **ratio of
discount rates**. (The example below adds a synthetic two-group factor for illustration.)

```{r emmeans, eval = has_tmb}
set.seed(2)
grp <- data.frame(
  id    = levels(jarvis2019$id),
  group = factor(sample(c("A", "B"), nlevels(jarvis2019$id), replace = TRUE))
)
jarvis_grp <- merge(jarvis2019, grp, by = "id")
fit_grp <- fit_dd_tmb(jarvis_grp, equation = "mazur", factors = "group")
get_dd_param_emms(fit_grp)                                  # EMM of k per group
get_dd_comparisons(fit_grp, contrast_by = "group", adjust = "holm")  # ratio-of-k contrasts
```

## Choosing a family

| | `family = "sltb"` (default) | `family = "gaussian"` |
|---|---|---|
| Respects $[0,1]$ | yes (incl. boundary mass at 0 and 1) | no (can predict outside) |
| Variance | shrinks near bounds, grows mid-range | constant |
| Boundary data | modeled | only via the bounded mean |
| Matches legacy `fit_dd()` | — | yes (NLS-style) |

Use `"sltb"` for inference on real IP data; use `"gaussian"` when you want continuity with a
classic least-squares discounting analysis or a quick baseline.

## How much to trust it

The SLT-beta density was checked before release: it integrates to one, has finite density at
both boundaries, reproduces the published moments, and recovers known discount rates from
simulated data (rank correlation 0.99 for Mazur, 0.95 for the exponential). On the Jarvis
(2019) data every boundary subject is fit successfully where ordinary beta regression cannot
be evaluated at all. The full Monte-Carlo study lives in the package's `dev/` directory.

## References

- Kim, M., Koffarnus, M. N., & Franck, C. T. (2024). *Thinking Inside the Bounds: Improved
  Error Distributions for Indifference Point Data Analysis and Simulation via Beta Regression
  Using Common Discounting Functions.*
- Kim, M., Kaplan, B. A., Koffarnus, M. N., & Franck, C. T. (2025). *Scale-Location-Truncated
  Beta Regression: Expanding Beta Regression to Accommodate 0 and 1.* arXiv:2509.13167.
- Jarvis, B. P., et al. (2019). Delay-discounting dataset (N = 126; 1 week – 8 years).
````

- [ ] **Step 2: Verify the vignette knits** (requires the `jarvis2019` fixture from Task S.4; API chunks only run if TMB is installed).

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'rmarkdown::render("vignettes/sltb-discounting.Rmd", quiet = TRUE)' && echo OK
```
Expected output: `OK` (an `.html` is produced in `vignettes/`). The `fit-tmb`/`emmeans` chunks evaluate when `has_tmb` is `TRUE`, and are skipped (shown but not run) otherwise. Remove the rendered `.html` afterward (`rm vignettes/sltb-discounting.html`) — it is a build artifact, not tracked.

- [ ] **Step 3: Apply Brent's manuscript voice on a final prose pass.** Invoke the `personalized-voice` skill (manuscript/scholarly voice) on the prose sections only (do not touch code chunks or the math). Re-knit to confirm it still builds.

- [ ] **Step 4: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add vignettes/sltb-discounting.Rmd && \
git commit -m "docs(vignette): ship sltb-discounting vignette on data(jarvis2019)

Converts the dev/ draft to vignettes/sltb-discounting.Rmd: conceptual + boundary
demonstration on the bundled Jarvis fixture, plus fit_dd_tmb/emmeans usage chunks
gated on TMB availability (CRAN-safe). Manuscript voice applied to prose.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task S.6: GitHub Actions R-CMD-check workflow with `build-essential`/`gfortran`

**Files:**
- Create: `/Users/brent/Dropbox/GIT/beezdiscounting/.github/workflows/R-CMD-check.yaml`
- Template (read-only): `/Users/brent/Dropbox/GIT/beezdemand/.github/workflows/R-CMD-check.yaml`

beezdiscounting has **no `.github/workflows/` directory** today. Port beezdemand's workflow verbatim, then add a Linux toolchain step so TMB compiles from source (TMB pulls C++/Fortran via RcppEigen). beezdemand's runner already had the toolchain implicitly; we add it explicitly to be safe across the matrix.

- [ ] **Step 1: Create the workflow** (port of `/Users/brent/Dropbox/GIT/beezdemand/.github/workflows/R-CMD-check.yaml` with one added Linux apt step before `setup-r-dependencies`):

```yaml
# Workflow derived from https://github.com/r-lib/actions/tree/v2/examples
# Need help debugging build failures? Start at https://github.com/r-lib/actions#where-to-find-help
on:
  push:
    branches: [main, master, develop, feat/*, fix/*]
  pull_request:
    branches: [main, master, develop]
  workflow_dispatch:

name: R-CMD-check

permissions: read-all

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}

    name: ${{ matrix.config.os }} (${{ matrix.config.r }})

    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: macos-latest,   r: 'release'}
          - {os: windows-latest, r: 'release'}
          - {os: ubuntu-latest,  r: 'devel', http-user-agent: 'release'}
          - {os: ubuntu-latest,  r: 'release'}
          - {os: ubuntu-latest,  r: 'oldrel-1'}

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
      OMP_NUM_THREADS: 1
      OPENBLAS_NUM_THREADS: 1
      MKL_NUM_THREADS: 1
      BLIS_NUM_THREADS: 1
      VECLIB_MAXIMUM_THREADS: 1
      R_DATATABLE_NUM_THREADS: 1

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
          http-user-agent: ${{ matrix.config.http-user-agent }}
          use-public-rspm: true

      - name: Work around ubuntu apt azure-cli repo 403
        if: runner.os == 'Linux'
        run: |
          set -euxo pipefail
          sudo rm -f /etc/apt/sources.list.d/azure-cli.list /etc/apt/sources.list.d/azure-cli.sources || true
          sudo apt-get update

      - name: Install C++/Fortran toolchain for TMB (Linux)
        if: runner.os == 'Linux'
        run: |
          set -euxo pipefail
          sudo apt-get install -y build-essential gfortran make

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check

      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
          build_args: 'c("--no-manual","--compact-vignettes=gs+qpdf")'
```

- [ ] **Step 2: Validate the YAML parses.**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'invisible(yaml::yaml.load_file(".github/workflows/R-CMD-check.yaml")); cat("YAML OK\n")'
```
Expected output: `YAML OK`.

- [ ] **Step 3: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add .github/workflows/R-CMD-check.yaml && \
git commit -m "ci: add R-CMD-check workflow with build-essential/gfortran for TMB

Ports beezdemand's R-CMD-check matrix (macOS/Windows/Linux x devel/release/
oldrel) and adds an apt step installing build-essential, gfortran, and make on
Linux so the TMB template compiles from source during the check.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task S.7: NEWS.md entry, DESCRIPTION → 0.4.0, `.Rbuildignore`, and document()/NAMESPACE regen

**Files:**
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/NEWS.md:1`
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/DESCRIPTION:5` (Version), plus `Imports`/`LinkingTo`/`Suggests`/`VignetteBuilder`/`SystemRequirements`
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/.Rbuildignore` (compiled artifacts, `dev`, `docs`)
- Modify: `/Users/brent/Dropbox/GIT/beezdiscounting/NAMESPACE` (regenerated by roxygen)

Note: the spec assigns the `DESCRIPTION` build-field edits (`LinkingTo: TMB, RcppEigen`; `Imports: TMB (>= 1.9.0), RcppEigen, emmeans`; `Suggests: knitr, rmarkdown`; `VignetteBuilder: knitr`; `SystemRequirements: GNU make`) to the TMB-engine phase's compile gate. If those edits already exist when this task runs, this task only **adds** `Version -> 0.4.0`, `Date`, the `Suggests: knitr, rmarkdown` (for the vignette) and confirms the fields are present. Do not duplicate fields — read `DESCRIPTION` first and add only what is missing.

- [ ] **Step 1: Bump the version and date, ensure vignette/build fields exist.** Edit `/Users/brent/Dropbox/GIT/beezdiscounting/DESCRIPTION`. Set:
  - `Version: 0.4.0`
  - `Date: 2026-06-06`
  - Confirm (add if missing) under `Imports:`: `TMB (>= 1.9.0)`, `RcppEigen`, `emmeans`.
  - Confirm (add if missing) under `Suggests:`: `knitr`, `rmarkdown` (testthat already present).
  - Confirm `LinkingTo: TMB, RcppEigen`, `VignetteBuilder: knitr`, `SystemRequirements: GNU make` exist.

Apply the version/date edit explicitly:
```
# old (DESCRIPTION lines 5-6):
#   Version: 0.3.2
#   Date: 2025-01-08
# new:
#   Version: 0.4.0
#   Date: 2026-06-06
```
Use the Edit tool to replace exactly those two lines. For the `Suggests:` block, append `knitr,` and `rmarkdown,` (alphabetical, before `testthat`) only if absent.

- [ ] **Step 2: Add the `.Rbuildignore` entries** for compiled objects, `dev/`, and `docs/`. Append to `/Users/brent/Dropbox/GIT/beezdiscounting/.Rbuildignore` (skip any line already present):

```
^src/.*\.o$
^src/.*\.so$
^src/.*\.dll$
^src/symbols\.rds$
^dev$
^docs$
```

- [ ] **Step 3: Write the NEWS.md entry** at the top of `/Users/brent/Dropbox/GIT/beezdiscounting/NEWS.md`. Insert a new `# beezdiscounting 0.4.0` block above the current development-version block:

```markdown
# beezdiscounting 0.4.0

### New Features

- **Mixed-effects discounting via TMB** (`fit_dd_tmb()`): fits the
  indifference-point (IP) family discounting model — Mazur hyperbolic or
  exponential mean with a subject random intercept on `log k` — under either the
  scale-location-truncated beta (`family = "sltb"`, default) or Gaussian
  (`family = "gaussian"`) observation family. Between-subject factors and
  continuous covariates enter the `log k` fixed-effect design.

- **SLT-beta error distribution**: assigns finite probability to indifference
  points at exactly 0 and 1, where ordinary beta regression is undefined. Means
  use an identity link on the discounting function; variance shrinks near the
  bounds and grows mid-range.

- **Estimated marginal means and contrasts**: `get_dd_param_emms()` returns the
  EMM of `k` per factor level (computed on the `log k` scale and
  back-transformed); `get_dd_comparisons()` returns pairwise or
  treatment-vs-control contrasts as ratios of discount rates, with multiplicity
  adjustment via any `stats::p.adjust` method.

- **broom + base S3 surface** on `beezdiscounting_tmb` objects: `tidy()`,
  `glance()` (backend `"TMB_mixed"`), `augment()`, `coef()`, `fixef()`,
  `ranef()`, `confint()`, `predict()`, `summary()`, `logLik()`, `AIC()`,
  `BIC()`, `nobs()`, `print()`.

- **`jarvis2019`** dataset: a de-duplicated Jarvis (2019) delay-discounting
  fixture (126 subject-records x 7 delays) with boundary indifference points,
  used by the new `sltb-discounting` vignette to demonstrate the SLT-beta
  family.

### Documentation

- New vignette `sltb-discounting`: why bounded error distributions matter for
  indifference points, the SLT-beta density, a boundary demonstration on the
  Jarvis (2019) data, and the mixed-effects workflow.

### Notes

- The data validator coerces percent/amount response scales to `[0, 1]` and
  clamps mild out-of-range values, **warning loudly** and naming the number of
  values coerced or clamped.

```

(Leave the existing `# beezdiscounting (development version)` and `# beezdiscounting 0.3.2` blocks below this new block.)

- [ ] **Step 4: Regenerate documentation and NAMESPACE.**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::document()'
```
Expected output: roxygen writes `man/*.Rd` (including `man/jarvis2019.Rd`, `man/fit_dd_tmb.Rd`, etc. produced by earlier phases) and updates `NAMESPACE` with the exports/imports from the new roxygen tags. No errors.

- [ ] **Step 5: Confirm the package builds and checks cleanly with the new fields.**

Run:
```
cd /Users/brent/Dropbox/GIT/beezdiscounting && Rscript -e 'devtools::check(document = FALSE, vignettes = TRUE, error_on = "warning")'
```
Expected output: `0 errors | 0 warnings | <n> notes` (TMB compiles; the vignette builds; version is 0.4.0). Investigate any error/warning before proceeding.

- [ ] **Step 6: Commit.**

```
cd /Users/brent/Dropbox/GIT/beezdiscounting && \
git add NEWS.md DESCRIPTION .Rbuildignore NAMESPACE man && \
git commit -m "chore(release): NEWS + DESCRIPTION 0.4.0 + document/NAMESPACE regen

Bump to 0.4.0; NEWS entry for the TMB mixed-effects discounting tier (fit_dd_tmb,
SLT-beta family, EMMs/contrasts, broom+base S3, jarvis2019 fixture, vignette);
.Rbuildignore compiled artifacts and dev/docs; regenerated man/ and NAMESPACE via
devtools::document().

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
