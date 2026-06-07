# Phase-1 Audit Findings — beezdiscounting mixed-effects discounting MVP (v0.4.0)

**Date:** 2026-06-07 · **Scope:** the IP-family mixed-effects discounting feature shipped on
`develop` (42 commits, unpushed). **Method:** three independent passes, cross-verified —
(1) a full manual code read of the entire feature; (2) authoritative `R CMD check --as-cran`
(`NOT_CRAN=true`, vignette + donttest run); (3) two adversarial reviewers — an 8-seam
multi-agent workflow with per-finding refutation (24 agents) and a Codex CLI read-only pass.
The workflow and Codex **converged** on the same defects; several were numerically reproduced.

---

## Verdict

**The statistical core is correct.** Both adversarial passes and the manual read agree:

- SLT-beta + Gaussian densities (`R/dd-density.R` ↔ `src/MixedDiscounting.h`) match the spec §3
  formula term-for-term; C++ `−nll` == R density to **2.3e-13**; the truncation normalizer `Z`
  is present and load-bearing; boundaries finite at `y=0`/`y=1`. `test-tmb-compile-gate.R`
  genuinely pins this tie-out for both equations × both families.
- C++ linear predictor `log_k_i = X.row(i)·β_k + σu·u_i`; subject-`k` reconstruction uses the
  subject's **design row** (factor-correct) and indexes **by name** (lexicographic alignment trap
  avoided); `predict()` looks up `u_i` by id name; phi floor applied on **both** optimizer paths.
- EMM / contrast **sign and scale** (log → log10 ÷ ln10, ratio = exp) are correct for full-rank,
  estimable grids. AIC/BIC df = `length(opt$par)` (β_k + σu + aux) — correct for a Laplace marginal.
- `jarvis2019` fixture is internally consistent and does **not** mis-group subjects.

The defects below are **contract / robustness / packaging** gaps, not math errors. None silently
corrupt a standard `id/x/y`, single-design, full-rank fit. They bite at the edges the per-task
reviews never exercised: non-default column names, non-PD Hessians, unbalanced interaction
designs, within-subject predictors, and CRAN packaging.

---

## P0 — real bugs that ship today (fix before release)

### B1. Non-default `y_var`/`x_var`/`id_var` break every response-path method
*Confirmed by both tracks; numerically reproduced.* `R/dd-tmb.R` `fit_dd_tmb` stores the user's
**original** column names in `param_info$x_var/y_var/id_var` (lines ~906–908), but `fit$data`
(`prepared$data`) is renamed to **canonical** `id/x/y` by the validator. So for any remapped fit —
e.g. `fit_dd_tmb(d, y_var="indiff", x_var="delay", id_var="subj")` — `predict(fit)` (newdata=NULL)
evaluates `object$data[["delay"]]` → `NULL`, and the subject path aborts on the missing `subj`
column. **`predict`, `fitted`, `residuals`, `augment` all fail.** The remap args are documented,
exported features, so this is a real broken-feature bug for anyone who uses them.
**Fix (simplest, complete):** store canonical names in `param_info` (`x_var="x"`, `y_var="y"`,
`id_var="id"`) — `fit$data` is always canonical and the by-name lookups then work; the examples
already pass canonical-named `newdata`. Add a regression test fitting with all three remapped.

### B2. `se_available` / `pdHess` gate is recorded but not enforced by SE-consuming methods
*Confirmed by both tracks (three loci).* The capstone fix set `se_available = sdr && hessian_pd`,
but downstream code ignores it:
1. `.dd_tmb_model_se()` (`R/dd-tmb-methods.R:54`) early-returns `model$se` **before** checking
   `se_available`, so `confint()` emits finite Wald CIs from a non-PD-Hessian fit with **no
   warning** (unlike `tidy()`/`summary()`, which attach one).
2. `get_dd_param_emms` / `.dd_compare_k` use `sdr$cov.fixed` whenever present, regardless of
   `se_available`, with no warning.
3. **Worst sub-point:** when `sdreport` fails, both EMM and contrast paths fall back to a
   **diagonal** vcov built from squared SEs — silently discarding off-diagonal covariance. A
   contrast variance `dxᵀ V dx` genuinely needs the off-diagonal terms, so contrast SEs/p-values
   become wrong (not just unavailable) on that path, with no warning.
**Fix:** gate all SE/CI/p-value output on `se_available`; NA-out or warn when false; on the
diagonal fallback, warn that contrast uncertainty is unreliable.

### B3. Rank-deficient / empty-cell interaction designs are accepted → meaningless EMMs
*Codex BLOCKING; workflow rated observation — real, narrow trigger.* `.dd_tmb_build_design`
sends `model.matrix()` straight to TMB with no rank check. With `factor_interaction=TRUE` and an
**empty cell** (unbalanced design, e.g. no rows with `A2 & B2`), the `A2:B2` column is all-zero
(aliased); TMB "estimates" an unidentified coefficient and `.dd_build_emm_ref_grid` averages over
the phantom cell in the full factorial crossing — producing a confident-looking but non-estimable
EMM. (The fit's Hessian is usually non-PD here, which *would* be flagged — but per B2 that gate
isn't enforced on the EMM path.) **Fix:** after design construction, `abort/warn if
qr(X)$rank < ncol(X)`, naming the aliased columns. One cheap check kills the whole class.

---

## P1 — contract guards (silent-surprise misuse; cheap to add)

### B4. Within-subject-varying predictor → arbitrary per-subject `k`; between-subject contract unenforced
*Both tracks (Codex major / workflow minor, confirmed).* The model fits on all rows (a valid
row-varying GLM), but `subject_pars$k`, `ranef()`, and `predict(type="parameters")` reconstruct
each subject's `k` from the **first** design row only (`R/dd-tmb.R:679`). If a declared
factor/covariate varies within `id`, the reported per-subject `k` is arbitrary, silently.
**Fix:** validate each declared predictor is constant within `id`; `cli_warn` (or `cli_abort` for
the documented between-subject MVP contract) naming offending ids/columns.

### B5. Declared `factors` not coerced; numeric "factor" silently fits as a continuous slope
*Codex.* `build_fixed_rhs` only drops single-level columns that are already `is.factor()`. A
numeric column passed as `factors=` is fit as a **continuous slope** (1 column, not level dummies)
— a silently different model — and char single-level factors aren't dropped. **Fix:** coerce
declared `factors` to factor (and drop/abort single-level) before fitting; require
`continuous_covariates` numeric + finite after complete-casing.

### B6. Canonical-name collision can overwrite `id`/`x`/`y`
*Codex (edge).* `.dd_validate_ip` doesn't reserve the canonical output names; `id_var="pid"` with
a factor literally named `"id"` overwrites the subject id via the extra-cols loop. **Fix:** treat
`id/x/y` as reserved; abort if an `extra_cols` name collides after remapping.

### B7. Blowup guard checks only `beta_k[1]` and only on the multi-start path
*Both tracks.* `.is_blowup` inspects the intercept alone; a degenerate non-reference level /
interaction / covariate slope can blow up observed `log k` while `beta_k[1]` stays sane, and
`multi_start=FALSE` has no guard at all. **Fix:** guard `max(abs(X %*% beta_k))` over observed
rows on **every** optimizer path. (The phi floor already blocks the main `k→∞`/`φ→0` sink, so this
is a secondary net.)

### B8. `predict()` raises an opaque base-R error when `newdata` omits a factor column entirely
*Completeness critic.* The unseen-level guard only runs when the column is present, so a `newdata`
that drops a needed factor/covariate falls through to an opaque `model.matrix` error. **Fix:**
check required columns up front; `cli_abort` cleanly.

---

## P2 — minor correctness / polish

- **B9. Percent auto-detect from a single stray value.** `max(y) > 1.5` divides the *whole* column
  by 100; one data-entry slip (`1.51`) silently rescales valid proportions. *Confirmed.* Consider a
  quantile/majority rule or an "ambiguous → error/ask" branch. (`R/dd-validate.R:146`)
- **B10. `.dd_normalize_re` accepts `log(k) ~ 1` / `I(k) ~ 1`** as if `k ~ 1` (`all.vars` strips
  wrappers). *Confirmed.* Tighten to require the bare `k` symbol. (`R/dd-re-utils.R:37`)
- **B11. `report_space = "log"` not publicly requestable** though the transform + registry support
  it (`tidy`/`summary`/`confint` enums omit it). *Partial.* Add the enum value or document the
  omission.

---

## CRAN / R-package readiness (Phase B) — from `R CMD check --as-cran`: 1 ERROR, 1 WARNING, 4 NOTEs

- **C1 — ERROR (blocking):** the 7 S3-method `\donttest` examples call the **internal**
  `.simulate_dd_ip_mixed()` → `could not find function` under `--run-donttest` (which CRAN runs).
  *(The handoff's "0 errors" was a plain check that skips donttest.)* **Fix options below — needs a
  decision.**
- **C2 — NOTE:** `Depends: R (>= 2.10)` understates the floor — native pipe `|>` ⇒ **R ≥ 4.2.0**,
  RDX3 `data/*.rda` ⇒ R ≥ 3.5.0. Bump to `R (>= 4.2.0)`. *Confirmed.*
- **C3 — NOTE:** `.Rbuildignore` omits `.claude/` and `personal_tests/` (both ship → hidden-dir +
  non-standard-top-level NOTEs). Add `^\.claude$` and `^personal_tests$`. *Confirmed.*
- **C4 — NOTE:** `NEWS.md` has an unparseable `# beezdiscounting (development version)` heading
  *below* the 0.4.0 section. Fold into 0.4.0 or retitle to a real version.
- **C5 — NOTE:** `URL https://brentkaplan.github.io/beezdiscounting/` 404s — expected until GitHub
  Pages is enabled (Phase F). Tracked, not a code fix.
- **C6 — WARNING:** install warning `-Wfixed-enum-extension` originates in R's own
  `R_ext/Boolean.h` under Apple clang 21 — **environmental**, not a package defect. Won't appear on
  CRAN's toolchain. Document in cran-comments.
- **C7 — INFO:** 13.6 Mb install (12.7 Mb TMB libs) — expected for a compiled package; note it.
- **Smaller:** `get_dd_param_emms`/`get_dd_comparisons` (+ `long_to_wide_mcq`/`wide_to_long_mcq`)
  lack `@examples`; DESCRIPTION carries both `Authors@R(cre)` and a manual `Maintainer:`; `inn()`
  documented but neither exported nor `\keyword{internal}`; `dd_ip` data doc has a vague `\format`.

---

## Test-coverage gaps (Phase A — backfill as regression tests)

All *confirmed*. These seams have no end-to-end test, which is why the contract bugs above slipped:

1. **factors + continuous covariate together** — never fit end-to-end.
2. **≥3 factors** — only `build_fixed_rhs` formula-level; never through fit/EMM/predict.
3. **`at=` restricting an omitted (marginalized) factor** — documented behavior, untested.
4. **`predict()` value on a non-reference factor cell** — only the unseen-level *error* is tested.
5. **`trt.vs.ctrl` × `contrast_by`** — never combined.
6. **unbalanced / empty interaction cell** — untested (this is the B3 trigger).
7. **`jarvis2019` fit end-to-end** (spec §4.8 boundary/degenerate regression) — only shape-tested.
8. Regression tests for the three *previously-fixed* seam bugs exist only partially (the 3rd-factor
   fix is asserted at formula level only).
9. Recovery asserts: prefer explicit `expect_lt(abs(est-truth)/truth, tol)` over
   `expect_equal(tolerance=)` for small targets (waldo's tolerance isn't relative there).

---

## The `jarvis2019` question (Brent, 2026-06-07: bundling not yet decided)

This gates several items above. If `jarvis2019` is **not** bundled:
- Examples must not use it (so the **C1 donttest fix should be self-contained** — export the
  simulator, or use a tiny inline `data.frame`).
- The Phase-E vignette's boundary demonstration needs a fallback (simulate boundary-heavy data, or
  ship a small purpose-built fixture).
- Test gap #7 (fit jarvis end-to-end) becomes "fit a simulated boundary-heavy dataset" instead.

If it **is** bundled, the above can lean on it. Either way, exporting the simulator as a public
`simulate_dd_ip()` makes examples/tests independent of that decision and is a listed deliverable.

---

## Recommended remediation sequencing (TDD, feature branch off `develop`, Codex-reviewed plan)

1. **CRAN-clean to a passing `R CMD check --as-cran`** (fast, unblocks everything):
   C1 (donttest), C2 (R floor), C3 (.Rbuildignore), C4 (NEWS), + the small doc items.
2. **P0 correctness** (one branch, TDD): B1 (canonical names), B2 (se gate, 3 loci),
   B3 (rank check) — each with a failing test first.
3. **P1 contract guards:** B4–B8.
4. **P2 polish:** B9–B11; backfill the 9 coverage gaps.
5. Re-run `R CMD check --as-cran` + full suite; Codex capstone re-review; then Phase B sign-off →
   Phase C (modeling fast-follows) per the handoff.

---

## Resolution log — branch `fix/phase1-cran-and-p0` (2026-06-07)

**RESOLVED (committed, verified):**
- **C1–C4 + T5(partial):** `R CMD check --as-cran` went from **1 ERROR, 1 WARNING, 4 NOTEs** →
  **1 WARNING, 1 NOTE**. ERROR gone (simulator exported as `simulate_dd_ip()`, examples valid
  under `--run-donttest`); the `.claude`/`personal_tests`/NEWS NOTEs gone; `Depends: R (>= 4.2.0)`.
  Remaining: the URL-404 NOTE (Pages — **C5/Phase F**) and the clang `-Wfixed-enum-extension`
  WARNING (**C6**, environmental, R's own header on Apple clang 21 — not a package defect).
- **B1** (canonical names): `param_info` now stores `id/x/y` (+ `user_vars`); remapped fits'
  post-fit methods work and match the canonical fit. Caller `newdata` must use canonical names —
  now **documented** + **clean `cli_abort`** when a required canonical column is absent.
- **B2** (se gate): enforced in `.dd_tmb_model_se` (gate before cached SEs), `confint` (warn +
  NA), and the shared `.dd_resolve_beta_vcov()` (NA + single warning for EMM/contrast; diagonal
  fallback only for the `se_available && cov.fixed=NULL` edge). Point estimates preserved.
- **B3** (rank guard): `.dd_tmb_build_design` aborts on `qr(X)$rank < ncol(X)`, naming aliased
  columns.

**Verification:** full local suite **667 pass / 0 fail / 0 warn / 0 skip**; `R CMD check --as-cran`
0 ERROR; Codex plan review **APPROVE-WITH-CHANGES** (incorporated) and Codex capstone re-review
**APPROVE-WITH-CHANGES** (incorporated: `newdata` doc + guards, warning wording).

## Resolution log — branch `fix/phase1-tranche2-p1-tests` (2026-06-07)

**RESOLVED (committed, verified):** the P1 contract guards + the audit's coverage-gap backfill.
- **B4** — warn (not abort) when a declared predictor varies within a subject
  (`.dd_check_between_subject`), naming the column(s) + affected-subject count.
- **B5** — coerce declared `factors` to factor (numeric "factor" no longer fits as a slope);
  require `continuous_covariates` numeric AND finite; reject predictors that name the
  id/delay/response column.
- **B6** — reserve canonical `id/x/y` in `.dd_validate_ip` (remapped-role collision).
- **B7** — generalize the blow-up guard to `max|X %*% beta_k|` (`.dd_logk_blowup`) on both the
  multi-start and single-start paths.
- **B8** — `predict()` errors cleanly when `newdata` omits a factor/covariate column.
- **Coverage gaps 1–7** backfilled as non-circular regression tests (factors+covariate together;
  3-factor additive AND interaction end-to-end; `at=` on a marginalized factor; predict on a
  non-reference cell; `trt.vs.ctrl`×`contrast_by` with per-by-cell p-adjust; unbalanced full-rank
  equal-weight averaging; boundary-heavy SLT-beta fit, jarvis-independent). Gap 8 (phi-floor on
  `multi_start=FALSE`) was already covered; gap 9 (relative-error asserts) is optional polish.

**Verification:** full local suite **711 pass / 0 fail**; `R CMD check --as-cran` unchanged
(1 env WARNING + 1 Pages NOTE); Codex plan + capstone reviews **APPROVE-WITH-CHANGES**
(incorporated). Merged `--ff-only` to `develop`.

**STILL DEFERRED:** P2 (B9 percent-detect robustness, B10 RE-LHS tighten, B11 `report_space="log"`);
the Phase-B CRAN-polish items (manual `Maintainer:`, `inn()`/`score_one_mcq27()` `@keywords
internal`, `dd_ip` `\format`, `@examples` on the mcq reshapers, spelling/urlchecker/lintr); gap-9
relative-error asserts. Then handoff phases B/C/D/E/F.
