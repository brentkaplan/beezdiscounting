# Plan — Phase-1 remediation: CRAN-clean + P0 correctness

**Date:** 2026-06-07 · **Repo:** beezdiscounting · **Base:** `develop` (42 commits unpushed)
**Branch:** `fix/phase1-cran-and-p0` (off `develop`; `--ff-only` merge back; delete after)
**Source:** `docs/audit-2026-06-07-phase1-findings.md`. **Scope (Brent-approved):** make
`R CMD check --as-cran` green + fix the three P0 correctness bugs (B1, B2, B3). P1/P2/test-gap
backfill are a later branch.
**Codex plan review (2026-06-07): APPROVE-WITH-CHANGES** — incorporated below: (R1) `confint()`
must NA-out uncertainty when `se_available=FALSE` (not warn-only); (R2) move the `se_available`
gate before the `model$se` early-return; (R3) when `!se_available`, EMM/comparison SE/CI/stat/p
are **NA** (the diagonal-fallback warning applies only to the `se_available=TRUE` but
`cov.fixed=NULL` edge); (R4) `data_all` is only *conditionally* canonical (the B6 name-collision
exception is out of scope here) — don't claim "always". Plus suggestions: dedup the inner EMM
warning in `.dd_compare_k`; add a tolerance note to the rank check; fold `score_one_mcq27()` into
T5 alongside `inn()`. **Decisions locked:** defer the `jarvis2019` bundling question
(nothing here depends on it); export the simulator as `simulate_dd_ip()`; examples use
`simulate_dd_ip()` or the already-bundled `dd_ip` (600×3 `id/x/y`, 125 exact-0/1 boundary points).

**Conventions:** TDD (BDD `describe()/it()`, testthat 3e), **test-first** for every behavioral
change; verification gate is real `R CMD check --as-cran` (`NOT_CRAN=true`), not `load_all`;
**targeted `git add <files>` only** (working tree has unrelated cruft — never `-A`/`.`); after any
roxygen change run `roxygen2::roxygenise()` (never hand-edit NAMESPACE — it drops `useDynLib`).
Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Manual
push only. TMB fits are slow (~min) — be patient; compiled tests need `skip_on_cran()` +
`skip_if_not_installed("TMB")`.

---

## Part 1 — CRAN-clean (target: 0 ERROR; only the expected Pages-404 NOTE + install-size INFO remain)

### T1 — Export the simulator as `simulate_dd_ip()` (fixes C1, the donttest ERROR)
- **Rename** `R/simulate-dd-mixed.R::.simulate_dd_ip_mixed` → `simulate_dd_ip` (same signature:
  `n_subjects, delays, log_k_pop, sigma_u, phi, sigma_e, family, equation, n_conditions, delta_k,
  seed`). Promote roxygen to a public doc: drop `@keywords internal`/`@noRd`, add `@export`,
  ensure `@param` for every arg, a `@return` describing the tibble, and a runnable `@examples`
  (small, fast: `simulate_dd_ip(n_subjects = 8, seed = 1)`).
- **Update all call sites** (mechanical) from `.simulate_dd_ip_mixed` → `simulate_dd_ip`:
  `R/dd-tmb-methods.R` (7 `\donttest` example blocks: predict/fitted/residuals/augment/confint/
  summary/print), `tests/testthat/{helper-dd-sim.R, test-simulate-dd-mixed.R, test-fit_dd_tmb.R,
  test-dd-tmb-methods.R, test-dd-comparisons.R}`.
- `roxygen2::roxygenise()` to regen NAMESPACE (`export(simulate_dd_ip)`) + `man/simulate_dd_ip.Rd`.
- **Test (first):** `test-simulate-dd-mixed.R` — add `it("is exported and returns the id/x/y(+condition) contract")`
  asserting `is.function(simulate_dd_ip)`, column names, row count, `y ∈ [0,1]`. The existing
  recovery tests (renamed) still pass. Confirm `"simulate_dd_ip"` ∈ `getNamespaceExports`.
- **Why rename vs wrapper:** one source of truth; tests then exercise the public API. No behavior
  change, so recovery tests are unaffected by the rename.

### T2 — DESCRIPTION R floor (C2)
- `Depends: R (>= 2.10)` → `R (>= 4.2.0)` (native pipe `|>` ⇒ 4.2.0; RDX3 data ⇒ 3.5.0; 4.2.0
  covers both). No test; verified by re-check (the build NOTE disappears).

### T3 — `.Rbuildignore` (C3) — stop shipping dev cruft
- Add: `^\.claude$`, `^personal_tests$`, and a global `\.DS_Store$` (there is a tracked
  `data/.DS_Store` and a top-level one). Also add `^\.Rhistory$`, `^\.Rprofile$` defensively.
- Verify via re-check: the hidden-dir + non-standard-top-level NOTEs disappear.

### T4 — NEWS.md (C4)
- The stray `# beezdiscounting (development version)` heading sits below `# beezdiscounting 0.4.0`.
  Since 0.4.0 is unreleased, **fold its `### Bug Fixes` content up into the 0.4.0 section** and
  delete the unparseable heading (cleanest; preserves the bug-fix notes). Re-check: the
  "Cannot extract version info" NOTE disappears.

### T5 — Small doc/metadata items (cheap CRAN hygiene)
- Add `@examples` to `get_dd_param_emms`, `get_dd_comparisons` (reuse a tiny `simulate_dd_ip()`
  fit under `\donttest`) and to `long_to_wide_mcq` / `wide_to_long_mcq` (plain `@examples`).
- DESCRIPTION: drop the manual `Maintainer:` line (CRAN derives it from `Authors@R` `cre`).
- `inn()`: add `@keywords internal` (+ `@noRd` if it should not get an Rd) — it is documented but
  neither exported nor internal-marked.
- `dd_ip` data doc: expand `@format` with a `\describe{}` of `id`, `x` (delay), `y` (proportion).
- `roxygenise()` after these. (These are "reduce NOTEs / reviewer-friendly"; if any proves
  fiddly under time pressure it can drop to the P2 branch — but all are one-liners.)

### T6 — Verify Part 1
- `R CMD build` + `R CMD check --as-cran --no-manual` with `NOT_CRAN=true`. **Gate:** 0 ERROR;
  remaining acceptable items are only the URL-404 (Pages, Phase F — document in `cran-comments.md`)
  and the 13.6 Mb install INFO (document). The clang `-Wfixed-enum-extension` WARNING is
  environmental (R's header on Apple clang 21) — note in cran-comments, not fixable in-package.

---

## Part 2 — P0 correctness (TDD: failing test → fix → green; no regressions)

### T7 — B1: canonical column names in `param_info`
- **Test (first), `test-dd-tmb-methods.R`:** fit two models on the same data — one default
  `id/x/y`, one remapped (`y_var="indiff", x_var="delay", id_var="subj"` on a renamed copy). Assert
  `predict(fit_remap)`, `fitted(fit_remap)`, `residuals(fit_remap)`, `augment(fit_remap)` all run
  and their `.fitted` equal the default fit's (same underlying data). Currently the remapped fit
  errors — RED.
- **Fix, `R/dd-tmb.R::fit_dd_tmb`:** set `param_info$x_var = "x"`, `y_var = "y"`, `id_var = "id"`
  (since `fit$data`/`data_all` are always canonicalized by `.dd_validate_ip`). Grep all
  `param_info$(x|y|id)_var` readers (`predict`, `.dd_tmb_predict_k`, `.dd_tmb_fitted_resid`) to
  confirm none need the original names. Document that `newdata` uses canonical `id/x/y`.
- *(Optional, only if trivial:)* also retain the user's original names under a separate
  `param_info$user_vars` for display — not required for the fix.

### T8 — B2: enforce the `se_available` / `pdHess` gate in SE-consuming methods (3 loci)
- **Tests (first), `test-dd-tmb-methods.R` + `test-dd-comparisons.R`:** construct/force a fit with
  `se_available = FALSE` (simplest: a fixture fit object with `se_available <- FALSE`,
  `hessian_pd <- FALSE`, and `sdr$cov.fixed` present-but-flagged). Assert:
  (a) `confint(fit)` returns **`NA` `conf.low/high`** (R1: not finite Wald + warning) and
  `cli_warn`s; (b) `get_dd_param_emms(fit)` / `get_dd_comparisons(fit)` emit **`NA`
  SE/CI/statistic/p** (R3) with a single `cli_warn`, not silent Wald CIs; (c) the
  `se_available=TRUE` but `cov.fixed=NULL` diagonal-fallback edge `cli_warn`s that contrast SEs
  drop off-diagonal covariance. RED on current code.
- **Fix:**
  - `R/dd-tmb-methods.R::.dd_tmb_model_se` — move the `se_available` check **before** the
    `model$se` early-return (R2): when `isFALSE(object$se_available)` return `NA` SEs, so
    `confint`/`tidy`/`summary` uniformly reflect the gate. `confint()` therefore yields `NA`
    intervals (R1). (`tidy`/`summary` already add the warning attr/note.)
  - `R/dd-comparisons.R::get_dd_param_emms` and `.dd_compare_k` (vcov at lines ~305 / ~547) —
    when `!isTRUE(fit$se_available)`, **do not** fall back to a diagonal vcov; set SE/CI/statistic/
    p to `NA` and `cli_warn` once (R3). The diagonal fallback (`diag(model$se^2)`) is retained
    **only** for the `se_available=TRUE` & `cov.fixed=NULL` case, with a warning that contrast
    variances `dxᵀ V dx` lose off-diagonal covariance. **Dedup:** `.dd_compare_k` calls
    `get_dd_param_emms` internally (line ~561) — suppress/share the warning so it fires once.
- **No-regression:** the existing happy-path EMM/contrast/confint tests (PD Hessian) must stay
  green — the gate only changes the non-PD / sdreport-failed paths.

### T9 — B3: rank-deficiency guard on the fixed-effect design
- **Test (first), new `test-dd-tmb-design-rank.R` (or in `test-fit_dd_tmb.R`):** build data with
  two factors + `factor_interaction = TRUE` and an **empty cell** (no rows with `A2 & B2`); assert
  `fit_dd_tmb(..., factor_interaction = TRUE)` aborts (or warns) naming the aliased column(s).
  Also a passing-case test: a full-rank 2-factor interaction design fits normally. RED on current.
- **Fix, `R/dd-tmb.R::.dd_tmb_build_design`:** after `X <- model.matrix(...)`, compute
  `r <- qr(X)$rank`; if `r < ncol(X)`, `cli_abort` (preferred — a rank-deficient design yields
  non-estimable coefficients) listing the aliased columns (those beyond the QR pivot rank). Keep
  the message actionable ("drop the interaction or collapse empty cells"). *Decision to confirm in
  Codex review:* abort vs warn-and-continue — abort is safer (the EMM path can't be trusted on a
  rank-deficient fit), but warn-and-fit-reference-only is an option. Default: **abort**.

### T10 — Verify Part 2 + whole
- Full suite via `NOT_CRAN=true Rscript -e 'pkgload::load_all(quiet=TRUE);
  dyn.load(TMB::dynlib("src/beezdiscounting")); testthat::test_dir("tests/testthat")'` → 0 fail.
- `R CMD check --as-cran` → 0 ERROR (Part-1 gate still holds).
- **Codex capstone re-review** of the diff (read-only) before merge.
- Merge `--ff-only` to `develop`; delete branch. **Do not push** (Brent's call). Update
  `docs/audit-2026-06-07-phase1-findings.md` to mark B1/B2/B3 + C1–C4 resolved.

---

## Risks / watch-items
- **B2 over-broadening:** NA-ing SEs must not break the many happy-path tests that assert finite
  SEs/CIs. Scope the gate strictly to `isFALSE(se_available)` / `sdr` NULL.
- **B1 newdata expectation:** after the fix, user `newdata` must use canonical `id/x/y`; the doc
  examples already do. Flag in roxygen so remap users aren't surprised.
- **T1 rename surface:** 7 example blocks + 5 test files; a missed call site fails the donttest
  check or a test — caught by T6/T10. Grep `\.simulate_dd_ip_mixed` must return zero after T1.
- **B3 abort vs warn:** confirm with Codex; abort is the safer default for a non-estimable design.
- Roxygen env: `roxygen2` is installed; `devtools`/`pkgdown`/`lintr` are not (install on demand).
