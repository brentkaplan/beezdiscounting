# Plan — Phase-1 remediation tranche 2: P1 contract guards + test-coverage backfill

**Date:** 2026-06-07 · **Repo:** beezdiscounting · **Base:** `develop` (tranche-1 merged, 46 ahead)
**Branch:** `fix/phase1-tranche2-p1-tests` (off `develop`; `--ff-only` back; delete after).
**Source:** `docs/audit-2026-06-07-phase1-findings.md` (P1 = B4–B8; the 9 coverage gaps).
**Scope (Brent-approved):** the P1 contract guards and the highest-value regression tests. P2
(B9–B11) is optional/stretch. **Same conventions as tranche 1:** TDD (BDD `describe/it`,
test-first), real `R CMD check --as-cran` gate, **targeted `git add` only**, `roxygenise()` after
roxygen changes, manual push, Codex-review this plan before coding, trailer
`Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

Context already in place from tranche 1: `param_info$factors` / `$continuous_covariates` are stored;
`fit$data` is canonical `id/x/y`; `simulate_dd_ip()` is exported.

**Codex plan review (2026-06-07): APPROVE-WITH-CHANGES — incorporated:**
- **B5:** also reject **non-finite** (`Inf`) covariates (`complete.cases` keeps `Inf`); and at
  `fit_dd_tmb` entry **reject `factors`/`continuous_covariates` that intersect the role vars**
  (`id_var`/`x_var`/`y_var`) — else `factors="x"` would coerce the delay column to a factor.
- **B4:** call after prepare **and after** B5's coercion, before design; use a **numeric tolerance**
  (warn "not exactly constant within id") so float noise doesn't false-warn; name the column(s) +
  affected-id count.
- **B7:** `.dd_logk_blowup(opt, X)` returns TRUE on length mismatch / non-finite `beta` / non-finite
  `eta` / `max(abs(eta)) > threshold`; thread into multi-start `.is_blowup` site **and** the
  single-start path after optimization.
- **Test #7:** use a **simulated boundary-heavy sltb fit** (fast, deterministic, jarvis-independent),
  NOT `dd_ip` (which has `y>1` clamp warnings + ~1.3 min runtime). **Test #8 DROPPED** — the phi
  floor on `multi_start=FALSE` is already asserted (`test-fit_dd_tmb.R:775`).

---

## Part 1 — P1 contract guards (each TDD: RED → fix → GREEN, no happy-path regression)

### B4 — between-subject contract: warn on a within-subject-varying predictor
A declared factor/covariate that varies within `id` still fits (valid row-varying GLM) but makes
`subject_pars$k` / `ranef()` / `predict(type="parameters")` reflect only the **first** design row.
- **Decision (confirm in Codex review):** **warn**, don't abort — the fit is valid; only the
  per-subject summary is approximate. `cli_warn` naming the offending column(s) + that subject-level
  `k` uses the first row. (Abort would reject a legitimate fixed-effects-only use of time-varying
  covariates.)
- **Fix:** a helper `.dd_check_between_subject(data, extra_cols, id_col = "id")` called in
  `fit_dd_tmb` after `.dd_tmb_prepare_data`; for each extra col,
  `tapply(col, id, \(v) length(unique(v)) > 1)` → warn if any subject varies.
- **Test:** build data where a covariate varies within one subject → `expect_warning(...,
  "within|between-subject")`; a clean between-subject covariate → `expect_no_warning`.

### B5 — declared `factors` are coerced; covariates required numeric
A numeric column passed as `factors=` silently fits as a **continuous slope** (1 col, not dummies);
char single-level factors aren't dropped (`build_fixed_rhs` only checks `is.factor`).
- **Fix (in `fit_dd_tmb`, on `prepared$data` before `.dd_tmb_build_design`):**
  `for (f in factors) prepared$data[[f]] <- as.factor(prepared$data[[f]])`; for each
  `continuous_covariates`, `cli_abort` if not numeric. Single-level drop in `build_fixed_rhs` then
  works for the coerced factors.
- **Test:** pass a numeric `grp` (values 1/2) as `factors="grp"`; assert the design has the factor
  **dummy** column (`grp2`), not a single continuous `grp` slope (recompute `ncol(X)` / colnames).
  A non-numeric covariate → `expect_error`.

### B6 — reserve the canonical `id`/`x`/`y` output names
`id_var="pid"` with a factor literally named `"id"` overwrites the canonical subject id in
`.dd_validate_ip` (`long[[col]] <- data[[col]]`).
- **Fix (in `.dd_validate_ip`, after computing `extra_cols`):**
  `collide <- intersect(extra_cols, c("id","x","y"))`; if non-empty, `cli_abort` ("rename the
  column; id/x/y are reserved canonical names").
- **Test:** `fit_dd_tmb(d, id_var="pid", factors="id")` → `expect_error(..., "reserved|canonical")`.

### B7 — generalize the optimizer blow-up guard to the full linear predictor
`.is_blowup` checks only `beta_k[1]` and only in the multi-start path; a degenerate non-reference
level / interaction / covariate slope can blow up observed `log k` while the intercept stays sane,
and `multi_start=FALSE` has no guard.
- **Fix:** a helper `.dd_logk_blowup(opt, X, max_abs = .dd_beta_k_abs_max)` computing
  `max(abs(X %*% beta_k))` over the observed design rows (the actual fitted `log k` range — avoids
  both the non-reference-level false-negative and the uncentered-covariate false-positive). Replace
  `.is_blowup` (pass `tmb_data$X`); add the same check (as a `cli_warn`) to the single-start path.
- **Test:** unit-test `.dd_logk_blowup` returns TRUE when a non-intercept `beta_k` pushes some row
  past the bound though `beta_k[1]` is sane; FALSE on a normal fit. (Avoid trying to force a real
  blow-up fit — test the predicate directly.)

### B8 — `predict()` clean error when `newdata` omits a needed factor/covariate column
The unseen-level guard runs only `if (f %in% names(newdata))`; a `newdata` omitting a factor/
covariate column entirely falls through to an opaque `model.matrix` error.
- **Fix (in `.dd_tmb_predict_k`, before `model.matrix`):** require all
  `c(pinfo$factors, pinfo$continuous_covariates)` present in `newdata`; `cli_abort` cleanly listing
  the missing ones. (Mirrors the tranche-1 `x`/`y` guards.)
- **Test:** fit with `factors="condition"`; `predict(fit, newdata=data.frame(x=c(7,365)))` (omits
  `condition`) → `expect_error(..., "condition|missing")`.

---

## Part 2 — regression test backfill (the audit's coverage gaps; some exercise Part-1 fixes)

New/expanded tests, all `skip_on_cran()` + `skip_if_not_installed("TMB")`, recomputing expectations
**non-circularly** (`beta %*% model.matrix`, never the function against itself):

1. **factors + continuous covariate together (end-to-end):** `fit_dd_tmb(sim, factors="condition",
   continuous_covariates="age")`; assert `subject_pars$k` and a per-cell EMM equal an independent
   recompute from `fit$formula_details$X` and `beta_k`.
2. **≥3 factors end-to-end:** fit with 3 factors (additive AND `factor_interaction=TRUE`); recompute
   one EMM cell non-circularly to confirm all three contribute (the B2-seam regression).
3. **`at=` restricting an OMITTED (marginalized) factor:** `get_dd_param_emms(fit, factors_in_emm="A",
   at=list(B="b1"))` equals the `B=b1`-only sub-average (hand-computed), not the average over all B.
4. **`predict()` value on a NON-reference factor cell:** `predict(level="population",
   newdata=<non-ref group>)` equals `mu(exp(beta_int+beta_grp), x)`, recomputed.
5. **`trt.vs.ctrl` × `contrast_by`:** assert row count `(nlev-1)*n_by`, each by-cell's reference is
   level 1, per-by-cell p-adjust.
6. **Unbalanced (but full-rank) cells:** a 2×2 with very different cell sizes fits and the EMM
   equal-weight averaging still matches the hand recompute (distinct from B3's empty-cell abort).
7. **Boundary-heavy end-to-end on `dd_ip`:** `fit_dd_tmb(dd_ip, family="sltb")` converges, `phi`
   stays ≥ the floor, population `k` finite/plausible — the jarvis-independent spec-§4.8 regression.
8. **phi-floor on `multi_start=FALSE`** (the third previously-fixed seam): assert a `multi_start=FALSE`
   sltb fit keeps `phi` ≥ `.dd_phi_min`.
9. *(stretch)* convert a couple of small-target recovery asserts to explicit
   `expect_lt(abs(est-truth)/truth, tol)`.

---

## Part 3 — verify
- Full suite via load_all recipe → 0 fail; `R CMD check --as-cran` → 0 ERROR (still only the
  Pages URL NOTE + env clang WARNING).
- `roxygenise()` if any roxygen changed (B5/B8 may add `@param`/notes).
- Codex capstone re-review of the diff; then `--ff-only` merge to `develop`, delete branch, mark
  B4–B8 + the gaps resolved in the audit doc. **No push** (Brent's call).

## Risks / decisions for Codex
- **B4 warn vs abort** (default: warn). **B7** predicate-only test (don't force a real blow-up).
- **B5 coercion site** (after prepare_data, before build_design) must not double-droplevels or
  desync `fit$data` from `X`. **B6** must not false-trigger for default `id/x/y` fits (extra_cols
  already excludes the role vars).
- Each Part-1 fix needs a happy-path no-regression assertion (the full suite covers most).
