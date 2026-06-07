# Design Spec — IP-family mixed-effects discounting (SLT-beta + Gaussian) via TMB

**Date:** 2026-06-06 · **Repo:** beezdiscounting · **Branch:** `feat/tmb-mixed-discounting`
**Scope of this push:** the **indifference-point (IP) family MVP** — 1-parameter discounting
(Mazur hyperbolic + exponential), random intercept on `log k`, between-subject factors,
the full S3 + emmeans/comparisons surface, simulation, and tests. The choice family
(Family 2), the 2-parameter tier (Green–Myerson / Rachlin), subject-random `φ`, and
parametric draws are **explicit fast-follows, out of scope here.**

This is "Plan A" of the 3-plan API initiative tracked in
`shinybeez/docs/plans/01-beezdiscounting-nlme.md`; it must establish the mixed-effects
discounting contract before the API (Plan B) is designed.

## 1. Methodological anchors (read before implementing)

- **Kim, Koffarnus & Franck (2024)** — *"Thinking Inside the Bounds"* — SLT-beta regression
  for hyperbolic-discounting indifference points, **per-subject** (no random effects).
  Mean uses an **identity link**: `μ_ij = 1/(1 + k_i·D_j)`, `k_i = exp(ψ_i)`. Reference R
  code at `/Users/brent/Dropbox/GIT/paper-code-kmg-slt_beta-47743946c512/R_code/`
  (`beta functions.R::scale_location_trun_beta_mazur`). Defines the IP simulation DGP.
- **Kim, Kaplan, Koffarnus & Franck (2025)** — arXiv 2509.13167 — generalizes SLTB and adds
  the **hierarchical** version (§3.3 *is* the mixed-effects discounting model: random
  intercepts on `ψ_i = log k_i` and `ln φ_i`). This MVP implements §3.3 with `φ` as a
  population parameter (subject-random `φ` deferred).
- **Port template:** beezdemand 0.3.0 TMB tier (`R/{param-space,param-registry,tmb-demand,
  tmb-methods,simulate-within-subject,tmb-parametric-draws}.R`, `src/{MixedDemand.h,
  beezdemand_TMBExports.cpp}`). The IP model is the **single-random-intercept special case**
  of that machinery, so most R scaffolding ports with simplification.

## 2. The model

Indifference proportion `y_ij ∈ [0,1]` for subject *i* at delay `D_j` (days).

- **Mean (identity link on a discounting function):**
  - Mazur (`equation = "mazur"`): `μ_ij = 1 / (1 + k_i·D_j)`
  - Exponential (`equation = "exponential"`): `μ_ij = exp(−k_i·D_j)`
  - μ is **guarded** to `[1e-6, 1−1e-6]` in the C++ template (exponential underflows to 0
    at long delays; verified necessary).
- **Subject discount rate:** `ψ_i = log(k_i) = Xβ + u_i`, `u_i ~ N(0, σ²)`.
  - `X` = fixed-effect design (between-subject factors + continuous covariates).
  - Single random intercept ⇒ no Cholesky/LKJ; `re_i = σ·u_i` (non-centered).
- **Observation family (`family =` dispatch):**
  - `family = "sltb"` (default): `y_ij ~ SLTBeta(μ_ij, φ)` on the closed `[0,1]`.
    `φ` = population precision, estimated as `log φ` (scalar). See §3 for the exact density.
  - `family = "gaussian"` (baseline): `y_ij ~ N(μ_ij, σ_e²)`, `σ_e` estimated as `log σ_e`.
    This is the 2024 paper's own comparator, the fast pipeline-bootstrap, and the cleanest
    IP↔NLS tie-out reference.

**Why the family is decoupled from emmeans:** the fixed/random effects sit on
`log k = Xβ` (linear in β). EMMs and contrasts are computed **for k** (the parameter of
scientific interest), on the `log k` scale, then back-transformed — exactly beezdemand's
"EMMs without emmeans" averaging-matrix trick. The family **only** changes the C++ density
term and the auxiliary parameter (`φ` vs `σ_e`); the entire emmeans/comparison surface is
identical across families.

## 3. The SLT-beta density (VERIFIED — implement exactly)

With `a = μφ`, `b = (1−μ)φ`, and **fixed constants `s = 1.0000001` (=1+1e-7), `l = 1e-8`**
(these are the *reference-code* values; the paper prose's `10^−8.5`/`10^−9` are a
prose-vs-code mismatch — use the code values):

```
log f(y | μ, φ) =  lgamma(a+b) − lgamma(a) − lgamma(b)            # = −log B(a,b)
                 + (a−1)·log(y/s + l) + (b−1)·log(1 − (y/s + l))   # scaled-beta kernel
                 − log(s)                                          # Jacobian
                 − log( pbeta(1/s + l, a, b) − pbeta(l, a, b) )    # truncation normalizer Z
```

- **Z is load-bearing**, NOT ≈1: `Z(a=0.05,b=1)=0.604`. The `−log(Z)` term must stay.
- **Boundaries are finite** at `y=0` (`y/s+l = 1e-8`) and `y=1` (`1−(1/s+l) = 9e-8`) — the
  intended feature; finite for all `a,b>0`.
- **TMB provides `pbeta`** (atomic, differentiable in `a,b`) — confirmed in TMB 1.9.19.
- **Moments** (for predict/augment): `E ≈ s(μ−l) ≈ μ`, `Var ≈ s²·μ(1−μ)/(φ+1)`. (Note the
  reference `var.y` drops the `s²`; we keep it / derive from SLT moments and document it.)
- **Mean–variance:** SLT/beta variance is **non-constant across delays** — a key advantage
  over Gaussian-on-y (which is homoscedastic and can predict outside `[0,1]`).

### Verification provenance (do not re-litigate; formalize as tests)
The density was verified three independent ways before this spec (scripts: `dev/sltb-
verification/`): (1) an independent mathematical re-derivation of the normalizer/Jacobian/
boundary behavior; (2) numerical property checks — normalization to 1e-13, kernel↔Z to
1e-14, boundary finiteness, beta-limit `f = dbeta/Z`, MC moments, exact match to the
reference body, MLE recovery (Mazur cor 0.986 / exponential 0.953 with multi-start);
(3) empirical fit to the Jarvis (2019) data (`Data_brent.csv`): **all 28 boundary subjects
fit with SLT where standard beta is non-finite; 126/126 converge; cor(log k)=0.96 vs NLS.**

## 4. Implementation landmines (from the reference audit — must address)

1. **Filter, never mask, out-of-range rows.** The reference `ll.temp * ifelse(IP∈[0,1],1,0)`
   is buggy: `0*NaN = NaN` poisons the whole subject. We never reach it because the
   validator coerces+clamps (§5), but the C++/R likelihood must also **drop NA rows up front.**
2. **Estimate `log k`, `log φ` (and `log σ_e`).** Unconstrained internal scale → positivity
   by construction; fixes the reference's unconstrained-optim NaN risk.
3. **Guard μ ∈ [1e-6, 1−1e-6]** in the template (exponential underflow; use `CppAD::CondExp*`,
   never `if` on Type-valued branches).
4. **Multi-start is essential**, especially for exponential (a fixed start gave cor=−0.26;
   multi-start gave 0.953). Port beezdemand's 3-set multi-start + data-driven starts.
5. **Simulator samples from SLT** (inverse-CDF via the truncated `pbeta`), not plain `rbeta`
   (the reference tutorial's plain-`rbeta` sim mismatches the fitted model at small shapes).
6. **Boundary asymmetry** (l added after ÷s ⇒ 1e-8 vs 9e-8 gaps). Reproduce the reference's
   exact scheme for tie-out; note a symmetric `[ε,1−ε]` mapping as a possible future refinement.
7. **`var.y` caveat:** the reference variance helper is the plain-beta variance (omits `s²`,
   ignores truncation). If we expose a variance helper, derive from SLT moments or document.
8. **Degenerate per-subject optimum (φ→0, k→∞).** Discovered during verification: for a
   boundary-heavy subject the SLT likelihood has a spurious global optimum where φ collapses
   to ~0 (all mass at the endpoints) with k→∞, scoring a *lower* nll than the sensible curve
   (real example: Jarvis subject 70, k 0.0038→414000, φ→0). Naive multi-start that keeps the
   global-min nll **selects** it. Mitigations: (a) the **random-effects prior on `log k`
   shrinks it away** — the principled cure and a core argument for the hierarchical model;
   (b) per-subject helpers (`.dd_tmb_compute_subject_pars`, `validate_subject_pars`) and the
   multi-start *selection* must guard against it — bound/penalize `log φ` away from −∞ (e.g.
   a weak prior or `φ ≥ φ_min`), and/or reject interior-implausible optima. Add a regression
   test for this subject/pattern.

## 5. Data contract & validator

- **IP family long format:** exactly `id, x, y` (matches beezdiscounting's existing
  `fit_dd`/`validate.R` convention; `x` = delay, `y` = indifference proportion). Expose
  `id_var/x_var/y_var` args to remap caller column names.
- **Response scale:** `y = IP/LL ∈ [0,1]`; **0 and 1 are valid** (SLTB's purpose).
- **Validator behavior (must WARN loudly and be DOCUMENTED — per Brent):**
  - percent inputs (detected as max > 1.5, or via an explicit `scale`/`ll` arg) divided by
    100 or by a supplied larger-later reward → `[0,1]`, **with a warning**.
  - mild out-of-range after coercion **clamped** (`>1 → 1`, `<0 → 0`) **with a warning that
    names how many values were clamped**. Never silent. Document in roxygen + a NEWS note.
- **Data caveat (Jarvis fixture):** `Data_brent.csv` `PID` is non-unique (78/126). If the
  package ships it as a fixture, dedupe to a unique subject id or document that rows are the
  subject unit; a *mixed* model treating duplicated PIDs as one subject would mis-group.

## 6. Files

**Create (R):** `R/dd-param-registry.R`, `R/dd-param-space.R`, `R/dd-tmb.R`
(`.dd_tmb_*` pipeline + `fit_dd_tmb()`), `R/dd-tmb-methods.R` (S3), `R/dd-comparisons.R`
(`get_dd_param_emms`/`get_dd_comparisons`), `R/simulate-dd-mixed.R`, `R/dd-validate.R`.
**Create (src):** `src/MixedDiscounting.h` (eqn_type × family), `src/beezdiscounting_TMBExports.cpp`
(`TMB_LIB_INIT R_init_beezdiscounting`), `src/Makevars`, `src/Makevars.win`.
**Create (tests):** `tests/testthat/helper-dd-sim.R`, `test-sltb-density.R` (formalize
`dev/sltb-verification`), `test-fit_dd_tmb.R`, `test-dd-tmb-methods.R`, `test-dd-comparisons.R`,
`test-dd-tmb-re-parser.R`, `test-simulate-dd-mixed.R`, `test-dd-nls-tieout.R`, plus backfill
`test-fit_dd.R` (the currently-untested NLS path).
**Create (vignette + data):** `vignettes/sltb-discounting.Rmd` — the user-facing
"what the package is doing and why" vignette (see §10); a **bundled, de-duplicated Jarvis
(2019) fixture** (e.g. `data/jarvis2019.rda` via `data-raw/`) so the boundary demonstration
runs on shippable data (the existing `dd_ip` fixture has few exact-0/1 boundaries; the Jarvis
data is the motivating case). Draft at `dev/sltb-verification/vignette-draft.Rmd`.
**Modify:** `DESCRIPTION` (`LinkingTo: TMB, RcppEigen`; `Imports: TMB (>= 1.9.0), RcppEigen,
emmeans`; `Suggests: knitr, rmarkdown`; `VignetteBuilder: knitr`; `SystemRequirements: GNU make`),
`NAMESPACE` (roxygen), `.Rbuildignore` (`^src/.*\.o$`, `^src/.*\.so$`, `^docs$`, `^dev$`),
`NEWS.md`, version → 0.4.0; add a CI workflow with `build-essential` (TMB compiles from source).

## 7. Build order (TDD — test named first)

1. `test-sltb-density.R` + `R/dd-param-space.R`/registry (pure R; port + SLT log-density in R) —
   density properties, delta-method round-trips. (No compile needed; fastest first slice.)
2. `R/dd-validate.R` + validator tests (coerce/clamp/warn).
3. `src/MixedDiscounting.h` + `TMBExports.cpp` + `Makevars*` + DESCRIPTION — compile gate.
   **Cross-check the C++ `-nll` against the R SLT log-density to 1e-8** (the standard TMB gate).
4. `R/dd-tmb.R` (`fit_dd_tmb()`, multi-start, sdreport + Hessian-PD gate) →
   `test-fit_dd_tmb.R` recovery (`exp(β̂₀)≈k_true`, tol 0.15) for both families × both equations.
5. `R/dd-tmb-methods.R` → `test-dd-tmb-methods.R` (broom contracts; `backend="TMB_mixed"`).
6. `R/dd-comparisons.R` → `test-dd-comparisons.R` (non-circular recompute; EMM↔contrast
   `log(10)` invariance).
7. `R/simulate-dd-mixed.R` → recovery on known truth (SLT + Gaussian DGP).
8. `test-dd-nls-tieout.R`: Gaussian family ties out to `nls`; SLTB `exp(β̂₀)` ≈ geometric-mean
   of two-stage k (NOT arithmetic mean — differs by `exp(σ̂²/2)`).

## 8. Validation gates

- `devtools::test()` green; `R CMD check` clean (TMB compiles).
- Density: C++ `-nll` == R SLT log-density to 1e-8; normalization, boundary finiteness, Z.
- Recovery within tolerance on simulated SLT *and* Gaussian data (both equations).
- IP↔NLS tie-out via **geometric mean**: `exp(β̂₀) ≈ median(two-stage k)` or
  `exp(β̂₀ + σ̂²/2) ≈ mean(two-stage k)`. Empirical anchor: Jarvis geometric-mean k = 0.00733.
- emmeans verified non-circularly (recompute from `coef %*% model.matrix`, never vs emmeans).

## 9. Risks

TMB compile from scratch (no `src/`/CI exists today; `TMB_LIB_INIT` must be
`R_init_beezdiscounting`). `pbeta` AD differentiability along the optimization path (verify
early). Identifiability of exponential at extreme k. `φ` (precision) identifiability with few
delays — kept population-level. **Degenerate φ→0/k→∞ optimum for boundary-heavy subjects**
(§4.8) — guard multi-start selection and subject-par computation; the RE prior is the cure.
Keeping the emmeans contract byte-compatible with beezdemand for the later unified API.

## 10. Vignette — `sltb-discounting.Rmd` (a first-class deliverable)

A user-facing vignette explaining *what the package is doing and why*. Outline:
1. **The problem** — discounting IPs are proportions in `[0,1]` that pile up at 0 (fully
   discounted) and 1 (no discounting); least-squares/Gaussian assumes unbounded,
   homoscedastic errors (can predict <0 or >1); standard beta regression is **undefined at
   exactly 0 or 1**.
2. **The SLT-beta solution** — scale-location-truncation gives finite mass at 0 and 1; the
   density (§3), with intuition (a microscopic stretch so the beta kernel is never evaluated
   at exactly 0/1) and the mean = discounting function (identity link), `φ` → variance that
   shrinks near the bounds and grows mid-range.
3. **Boundary demonstration on real data** (the bundled Jarvis fixture) — show the boundary
   IPs, show standard beta is non-finite on those subjects while SLT fits them, compare SLT-k
   to NLS-k. (This is the runnable "gold" — works today.)
4. **The mixed-effects model** — random intercept on `log k`, between-subject factors,
   `fit_dd_tmb()`, parameter interpretation, `get_dd_param_emms()`/`get_dd_comparisons()`.
5. **Choosing a family** — SLT-beta (default) vs Gaussian (baseline): when/why.
6. **Confidence** — a short recovery + boundary-advantage summary (from §3 verification).

Authoring: the conceptual + boundary-demo spine (1–3, 6) is drafted now and runs today; the
API-usage sections (4–5) are written as `fit_dd_tmb()`/emmeans land (flip chunks to
`eval=TRUE`). Apply Brent's **manuscript/scholarly voice** (personalized-voice skill) on the
final pass. Keep the heavy Monte-Carlo validation in `dev/`/a pkgdown article, not the CRAN
vignette (keep the shipped vignette lean).

## 11. Out of scope (fast-follows)

2-parameter tier (Green–Myerson, Rachlin via `eqn_type`; `s` on `log s`, started at 1, no hard
`s≤1` cap); choice family (Family 2: structural + descriptive, binomial GLMM); subject-random
`φ` (→ the real pdDiag/pdSymm); parametric draws; Ebert–Prelec; probability discounting.
