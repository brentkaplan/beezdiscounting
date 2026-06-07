# Plan — Phase B: CRAN-readiness polish

**Date:** 2026-06-07 · **Repo:** beezdiscounting · **Base:** `develop` (49 ahead, P0+P1 merged)
**Branch:** `chore/phase-b-cran-readiness`. **Type:** resubmission (CRAN history 0.1.0→0.3.2).
**Conventions:** as tranches 1–2 (TDD for behavior, targeted `git add`, `roxygenise()` after roxygen,
real `R CMD check --as-cran` gate, Codex-review this plan, manual push, trailer). Diagnostics already
run: `spelling`, `urlchecker`, `lintr` (installed to the user lib; run them from OUTSIDE the package
renv by passing the absolute path, e.g. `spelling::spell_check_package("/Users/.../beezdiscounting")`).

The statistical core + Title/Description/`[cph]`/`@return` coverage are already CRAN-compliant. This
phase is doc hygiene + a lint config + the P2 audit items.

**Codex plan review (2026-06-07): APPROVE-WITH-CHANGES — incorporated:**
- **B9:** `median(y)>1.5` is wrong (`c(100,0,0,0,0)` → median 0 → would abort genuine percent). Rule:
  on the auto branch, if `any(y>1.5)`, auto-percent **only when the fraction of POSITIVE values
  exceeding 1.5 is ≥ 0.5 AND `max(y) ≤ 100`**; otherwise **abort** ("set `response_scale`"). Update
  the stale `R/dd-validate.R:17` roxygen that promises `max(y)>1.5`.
- **B10:** `is.name()` alone is incomplete (`phi ~ 1` is a name). Use
  `lhs <- random_effects[[2]]; !is.symbol(lhs) || !identical(lhs, quote(k))` (accepts `k~1` and
  `` `k`~1 ``; rejects `log(k)~1`, `I(k)~1`, `phi~1`, `k+phi~1`). Add **direct parser tests** in
  `test-dd-tmb-re-parser.R`, not only a slow `fit_dd_tmb` integration test.
- **B11:** also add `"log"` to the roxygen `@param report_space` text and fix the
  `test-dd-tmb-methods.R:390` "all three report_space values" test.
- **README (NEW required):** README says 0.3.2 in several places, points users to
  `citation("beezdemand")`, and references `hex-beezdiscounting.png` which `.Rbuildignore` excludes
  (relative link to a build-missing file → CRAN file-URI flag). Edit **README.Rmd** (version→0.4.0,
  fix citation, fix/absolute-URL the hex image), then `devtools::build_readme()`.
- **Lint:** `subqids` (`R/mcq.R:273`) is genuinely dead → remove; `safe_fit_model`/`x_char` are NOT
  dead (keep); the cli-glue items (`not_in`/`unknown_ids`/`aliased`/`n_affected`) are false positives.
  Disabling `object_usage_linter` is a stated tradeoff (glue interpolation isn't parsed).
- **Suggested:** add the SLT-beta method refs (Kim 2024; Kim/Kaplan 2025, arXiv:2509.13167) to the
  DESCRIPTION Description (not just cran-comments).

---

## Part A — CRAN doc/metadata hygiene (mechanical; verified by R CMD check)

1. **Typos** (real, from spelling): `defualt`→`default` (`R/utils.R:50` roxygen), `ouutput`→`output`
   (`R/mcq.R:323` roxygen). roxygenise after.
2. **`inst/WORDLIST`**: after fixing the typos, `spelling::update_wordlist()` to capture the
   legitimate technical terms/author names (beezdemand, emmeans, EMM/EMMs, SLT/sltb, TMB, Mazur,
   Rachlin, Myerson, nlme, nlminb, sdreport, pbeta, Qualtrics column names, Koffarnus/Bickel/Jarvis/
   … author names, LaTeX `infty`/`tfrac`/`mathrm`/`underbrace`, arg names). Add `Language: en-US` to
   DESCRIPTION (spelling requested it).
3. **`@examples`** on the two exports missing them: `long_to_wide_mcq`, `wide_to_long_mcq` (use the
   bundled `mcq27` data or a tiny inline frame; plain `@examples`, fast).
4. **`inn()` / `score_one_mcq27()`** (documented-but-unexported): add `@keywords internal` (keep the
   page, de-index) — they are user-relevant helpers, so prefer `@keywords internal` over `@noRd`.
   Confirm both have `@return`.
5. **`dd_ip` `\format`**: expand the vague one-liner to a `\describe{}` of `id` (subject), `x`
   (delay), `y` (indifference proportion).
6. **Drop the manual `Maintainer:`** line in DESCRIPTION (auto-derived from `Authors@R` `cre`).
7. **`cran-comments.md`**: rewrite for 0.4.0 — test environments; R CMD check result (1 env
   WARNING + 1 NOTE); document (a) the ~12.7 Mb install size (compiled TMB libs, acceptable for a
   compiled pkg); (b) the `-Wfixed-enum-extension` install WARNING is from R's own header under Apple
   clang and not a package defect; (c) the `brentkaplan.github.io/beezdiscounting/` URL resolves once
   the pkgdown site deploys with this release (aspirational); (d) the doi.org 403s are valid DOIs that
   block automated HEAD requests; (e) method references for the new SLT-beta feature (Kim/Koffarnus/
   Franck 2024; Kim/Kaplan/Koffarnus/Franck 2025, arXiv:2509.13167) — add a one-line note or include
   in DESCRIPTION (decide in review).
8. **URLs:** no code change — github.io is aspirational (Phase F deploys it); doi 403s are valid.
   `urlchecker::url_check()` re-run after Pages is live (Phase F), not now.

## Part B — lint config + substantive lints

9. **Add `.lintr`** encoding the package's actual style (mirrors shinybeez's philosophy):
   `linters: linters_with_defaults(line_length_linter(120), object_length_linter(40),
   return_linter = NULL, commented_code_linter = NULL, object_usage_linter = NULL)` — the 602
   default-strict lints are 80-char/snake_case/implicit-return style the established package never
   used; `object_usage` is mostly FALSE positives here (cli-glue interpolation like `{.val {aliased}}`,
   and lintr not resolving the package's own internal fns). With the config, the genuine signal is
   small.
10. **Genuine dead-variable lints** (pre-existing, optional): `safe_fit_model` (`fitting.R:93`),
    `x_char` (`fivetrial.R:2,139`), `subqids` (`mcq.R:273`) — inspect; remove only if truly dead and
    low-risk. Do NOT touch the cli-glue false positives (`not_in`, `unknown_ids`, `aliased`,
    `n_affected`).

## Part C — P2 audit items (TDD: failing test → fix → green)

11. **B10 — tighten the RE LHS** (`R/dd-re-utils.R:38`): require the LHS to be the BARE symbol `k`
    (`is.name(random_effects[[2]])`), so `log(k) ~ 1` / `I(k) ~ 1` are rejected with the existing
    fast-follow message instead of silently treated as `k ~ 1` (`all.vars` strips wrappers).
    *Test:* `expect_error(fit_dd_tmb(..., random_effects = log(k) ~ 1), "k")`.
12. **B11 — expose `report_space = "log"`** in `tidy()` (`:615`) and `summary()` (`:842`) enums (the
    delta-method transform + registry already support it; `confint`'s `"internal"` is already log for
    `beta_k`, leave it). *Test:* `tidy(fit, report_space="log")` returns `beta_k` on the log scale
    (== internal estimates) with `estimate_scale == "log"`.
13. **B9 — robust percent auto-detect** (`R/dd-validate.R:180`): replace the bare `max(y) > 1.5`
    trigger. Rule: on the `proportion` (auto) branch, divide by 100 only when the data are clearly
    percent (`stats::median(y, na.rm=TRUE) > 1.5`); if `max(y) > 1.5` but the median is not (a few
    out-of-range outliers among proportions), **abort** asking the caller to set `response_scale`
    explicitly — never silently rescale a whole column from a stray value. *Tests:* genuine percent
    (0–100) still auto-divides + warns; one stray 1.51 among proportions now aborts (not silent
    ÷100); confirm the existing percent-detection tests still pass (they use 0–100 data, median>1.5).

## Part D — verify
- `roxygenise()`; full suite (load_all recipe) 0 fail; `R CMD check --as-cran` → still 1 WARNING
  (env clang) + 1 NOTE (aspirational URL); `spell_check_package` clean (WORDLIST); `lint_package`
  signal-only. Codex capstone of the diff; `--ff-only` merge to develop; mark resolved. No push.

## Decisions for Codex review
- **B9 abort-on-ambiguous** vs warn-and-clamp (default: abort — safest; never silently rescale).
- **B11** scope (tidy+summary only; confint already log via "internal").
- **`.lintr`** disabling `object_usage_linter` (mostly false positives here) vs keeping it + a NEWS
  of accepted noise — default disable, since the cli-glue/internal-fn false positives dominate.
- **cran-comments method references** — add the SLT-beta refs to DESCRIPTION Description, or only to
  cran-comments? (Description is already long; default: cran-comments note + leave Description.)
- Confirm dropping the manual `Maintainer:` is safe (Authors@R `cre` present) and that the
  established Title/Description need no change.
