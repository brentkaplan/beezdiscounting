# Design Spec — 2-parameter discounting tier (Green-Myerson + Rachlin)

**Date:** 2026-06-07 · **Repo:** beezdiscounting · **Branch:** `feat/2param-discounting-tier`
**Phase:** the first Phase-C modeling fast-follow (handoff §C.1; SLT spec §11 "out of scope").
**Scope:** add the two standard 2-parameter hyperboloid discounting equations (Green-Myerson and
Rachlin) to the existing TMB mixed-effects machinery, with the nonlinearity exponent `s` as a single
**population** parameter. The choice family (Family 2), subject-random `s`/`φ`, and parametric draws
remain separate later slices.

## 1. Methodological anchors

- **Mazur (1987)** hyperbola, already implemented: `μ = 1 / (1 + kD)`.
- **Green & Myerson (2004) / Myerson & Green (1995)** hyperboloid: `μ = (1 + kD)^(−s)`. The exponent
  `s` scales the *whole denominator*; `s` captures the rate of decline beyond the simple hyperbola.
  Reduces to Mazur at `s = 1`.
- **Rachlin (2006)** hyperboloid: `μ = 1 / (1 + k·D^s)`. The exponent `s` sits on *delay*. Reduces to
  Mazur at `s = 1`. (McKerchar et al., 2009, compare these forms on the same data.)
- Both 2-parameter forms are nested over Mazur: testing whether `log s` departs from 0 tests curvature
  beyond the simple hyperbola.

## 2. The model (what changes vs the 1-parameter MVP)

The subject discount rate, random effect, fixed-effect design, and the two observation families
(`sltb` default, `gaussian`) are **unchanged**. The mean function gains a second shape parameter:

- Mean (`equation =` dispatch, identity link):
  - `mazur` (`eqn_type 0`): `μ = 1 / (1 + k_i·D)`  *(unchanged)*
  - `exponential` (`eqn_type 1`): `μ = exp(−k_i·D)`  *(unchanged)*
  - `green-myerson` (`eqn_type 2`): `μ = (1 + k_i·D)^(−s)`
  - `rachlin` (`eqn_type 3`): `μ = 1 / (1 + k_i·D^s)`
  - `μ` stays guarded to `[1e-6, 1−1e-6]` (the GM/Rachlin forms can also approach the bounds).
- **`s` is a single population scalar** (NOT in the `log k = Xβ` design, NOT a random effect), estimated
  on the log scale: `s = exp(log_s)`, started at `log_s = 0` (`s = 1`). **No hard `s ≤ 1` cap** — only
  soft/wide optimizer bounds (e.g. `log_s ∈ [log(0.05), log(20)]` as a numerical guard, not a
  constraint of scientific interest). Random `s` is poorly identified and is a later slice.
- `k`, the random intercept on `log k`, the between-subject factors/covariates, and the family-specific
  auxiliary (`φ` / `σ_e`) all behave exactly as in the 1-parameter tier. The family only changes the
  error term; `s` only changes the mean. The two are orthogonal.

## 3. Parameterization & the `s` naming decision (locked with Brent)

- The fitted exponent takes the canonical field name **`s`** (what users, `emmeans`-style output, and
  display labels expect for the hyperboloid). The param registry's `s` entry is **repurposed** from the
  SLT density's scale constant to this discounting exponent (`canonical = "s"`, `constraint = "s > 0"`,
  `valid_scales = c("natural","log","log10")`, `default_scale = "natural"`).
- The SLT-beta density's scale micro-constant (`1.0000001`) was never a fitted parameter; it stays a
  **fixed, unnamed C++ constant** (`Type s_slt = 1.0000001;`) with **no registry entry**. (Rename the
  C++ local from `s` to `s_slt` inside the SLT branch so it does not shadow the discounting `s`.)
- **Also rename the R-side SLT-scale locals to `s_slt`** (Codex suggestion 1): the `s` argument/local in
  `R/simulate-dd-mixed.R` (`s <- 1.0000001`) and the `s`/`@param s` in `R/dd-density.R::.dd_slt_logpdf`
  refer to the SLT scale constant; rename to `s_slt` so that `s` unambiguously means the discounting
  exponent package-wide. These are internal locals/args (no API impact).

## 4. Implementation

### 4.1 C++ template `src/MixedDiscounting.h`
- **First** rename the SLT-beta branch's local scale constant `s` → `s_slt` (the `1.0000001`), freeing
  the name `s` for the discounting exponent (Codex 8: this is required so `ADREPORT(s)` reports `s`).
- Add `PARAMETER(log_s);` and `Type s = exp(log_s);` near the other scalar transforms.
- Extend the mean switch (currently `if (eqn_type == 0) {...} else {...}`) to a 4-way `if` on the
  **integer** `eqn_type` (it is a `DATA_INTEGER`, so a plain C++ `if` is allowed — only `Type`-valued
  conditions need `CppAD::CondExp`):
  - `0` (mazur): `mu_raw = 1/(1 + k_i*x(i))`
  - `1` (exponential): `mu_raw = exp(-k_i*x(i))`
  - `2` (green-myerson): `mu_raw = pow(1 + k_i*x(i), -s)`
  - `3` (rachlin): **explicit `x == 0` guard** (Codex 1 — `pow(0, s)`'s value is 0 but its derivative
    w.r.t. `s` involves `log(0)`):
    ```cpp
    Type x_safe  = CppAD::CondExpGt(x(i), Type(0.0), x(i), Type(1.0)); // never feed 0 to pow
    Type rachlin = Type(1.0) / (Type(1.0) + k_i * pow(x_safe, s));
    mu_raw       = CppAD::CondExpGt(x(i), Type(0.0), rachlin, Type(1.0)); // x=0 -> mu=1
    ```
    Do NOT use `max(x, tiny)` (with `s` near 0, `tiny^s` is not tiny). The validator permits `x = 0`
    (`R/dd-validate.R` rejects only `x < 0`), so this branch is reachable.
- `ADREPORT(s)` alongside the existing ADREPORTs. (For 1-parameter fits `log_s` is `map`-fixed in §4.2,
  so `s = 1` is a constant with no usable SE; the S3 layer surfaces `s` only when `has_s`.)
- `1 + k*x > 0` always and `1 + k*x_safe^s > 0`, so the `pow` bases are strictly positive and AD-safe.

### 4.2 R fit pipeline `R/dd-tmb.R`
- `fit_dd_tmb(equation = c("mazur","exponential","green-myerson","rachlin"))`.
- `.dd_tmb_build_tmb_data`: extend the `eqn_type` switch with `"green-myerson" = 2L, "rachlin" = 3L`.
- `.dd_tmb_default_starts`: add `log_s = 0` (s = 1 start). The data-driven `beta_k[1]` start keeps the
  Mazur inversion (a reasonable `k` start for all four equations).
- **Coefficient scale contract (Codex 2):** `log_s` is kept under its optimizer name `log_s` in
  `coefficients`/`se` (the log scale), exactly mirroring how `log_aux` is kept as `log_phi`/`log_sigma_e`
  — it is NOT stored as natural `s`. The natural-scale `s` and its CI come from the param-space transform
  in `tidy()`/`summary()`/`confint()` (display term `s`, `estimate_scale = "log"`; §4.3/§4.5). Predict/
  fitted read `exp(coef[["log_s"]])`. This avoids a double delta-transform.
- **The `map` (Codex 4):** a new helper `.dd_tmb_build_map(has_s)` returns `list(log_s = factor(NA))`
  when `!has_s` (1-parameter: hold `log_s` fixed, never estimated/unidentified) and `NULL` when `has_s`.
  Thread `map` through **EVERY** `TMB::MakeADFun(...)` call: the single-start path, each multi-start
  build, AND the two test files that build `MakeADFun` directly (`test-tmb-compile-gate.R`,
  `test-fit_dd_tmb.R`) once `log_s` is in the parameter list. A mapped scalar coexists with `random = "u"`
  and `sdreport`.
- `.dd_tmb_extract_estimates(obj, opt, n_subjects, family, has_s, verbose)` (Codex 5): pass `has_s` so it
  surfaces `log_s` in `coefficients`/`se` only for 2-parameter fits (for 1-parameter fits `log_s` is
  mapped, so it is absent from `opt$par`; treat the presence of a free `log_s` as authoritative and assert
  it matches `has_s`). Pull its SE from `summary(sdr, "fixed")` like the other scalars.
- **`start_values$s` alias (Codex suggestion 3):** before merging user `start_values`, if the user passes
  `s` (the public name), convert to `log_s = log(s)`; otherwise `start_values$s` is silently dropped (the
  override only honors names already in `default_starts`).
- Multi-start: the 3 sets perturb `beta_k[1]`/`log_sigma_u`; for 2-parameter equations also perturb
  `log_s` (one set `log(0.7)`, one `log(1.4)`), since `s` can be the poorly-identified axis. The blow-up
  guard (`.dd_logk_blowup`) is unchanged.
- Soft `log_s` optimizer bounds (e.g. `[log(0.05), log(20)]`) added to `tmb_control$lower/$upper` for
  2-parameter equations only via `.expand_bounds`; user values win (mirrors `.dd_apply_phi_floor`).
- `param_info` records `equation` (already) and `has_s = equation %in% c("green-myerson","rachlin")`.

### 4.3 S3 surface `R/dd-tmb-methods.R`
- `.dd_discount_mu(k, x, equation, s = 1)`: extend to GM/Rachlin (used by `predict`/`fitted`/`augment`).
  Default `s = 1` keeps 1-parameter callers unchanged; methods pass `s = exp(coef[["log_s"]])` for
  2-parameter fits (with the same Rachlin `x == 0 -> mu = 1` guard as the C++).
- **`.dd_tmb_build_term_names()` (Codex 3):** map the raw coefficient name `log_s` to the display term
  `s`. The param-space transformer keys on `^s($|_)` and transforms a term named `s` but NOT `log_s`, so
  this display rename is what makes `tidy(report_space=)` / `confint(report_space="natural")` exponentiate
  it correctly.
- **`tidy()` (Codex 6):** currently emits only `beta_k` rows as `component = "fixed"`. Add the `log_s`
  row for 2-parameter fits as a new `component = "shape"` row (display `s`, `estimate_scale = "log"`,
  Wald `statistic`/`p.value` on the log scale), then run it through `.dd_transform_coef_table` so
  `report_space` back-transforms estimate+SE. This stays within the exact 8-column broom contract (a new
  row, not a new column). **`s` is NOT a variance component** — those rows intentionally carry `NA`
  SE/stat/p; `s` has a real Wald SE.
- **`summary()` (Codex 6):** currently classifies every non-`beta_k` coefficient as `variance` and drops
  it from the coefficient table. Add explicit `log_s` handling so `s` appears in the fixed/shape
  coefficient table (transformed by `report_space`), not silently dropped; update the print header/notes.
- **`.dd_tmb_model_se()` (Codex 7):** the SE-name translation currently only maps `log_aux`; ensure the
  `match` covers `log_s` too.
- **`confint()` (Codex 7):** currently only `beta_k` rows are exponentiated under `report_space =
  "natural"`. Include the `log_s` row in the natural-scale back-transform (`s = exp(log_s)`) so its
  interval is on the natural `s` scale.
- `coef()` returns the flat `coefficients` vector, now including `log_s` for 2-parameter fits (the
  optimizer-scale escape hatch, unchanged in spirit). `glance()` is unchanged.
- **emmeans/contrasts are unchanged: `k` only.** No per-group `s` in this slice, so
  `get_dd_param_emms`/`get_dd_comparisons` are untouched.

### 4.4 Simulation `R/simulate-dd-mixed.R`
- `simulate_dd_ip(..., equation = c("mazur","exponential","green-myerson","rachlin"), s = 1)`: compute
  `mu` with the chosen equation and `s`; the SLT/Gaussian draw is unchanged. Gives a known-`(k,s)` truth
  for recovery tests.

### 4.5 Registry/param-space `R/dd-param-registry.R`, `R/dd-param-space.R`
- Repurpose the `s` registry entry (§3). The param-space transforms/display already key on `^s($|_)` and
  handle `s` generically, so `tidy(report_space=)` on `s` works once `s` is a reported coefficient.

## 5. Build order (TDD — test named first)

1. **R reference mean + density extension** + `test-sltb-density.R` additions: extend the R `.dd_*` mean
   helper to GM/Rachlin; assert each reduces to Mazur at `s = 1` (pure R, no compile).
2. **C++ template** + `test-tmb-compile-gate.R`: add GM/Rachlin gate cases (single-subject, intercept-
   only) asserting **C++ −nll == R density to 1e-8** for both new equations × both families, at `s ≠ 1`.
   The compile gate is the load-bearing correctness check (it pins the `pow` arithmetic, incl. a Rachlin
   case with a `x = 0` row). The gate's `.gate_setup` parameter list must include `log_s`, and its direct
   `MakeADFun` must `map` `log_s` for the 1-parameter cases (leave it free for GM/Rachlin) — see Codex 4.
3. **Fit pipeline + map** + `test-fit_dd_tmb.R`: (a) the `map` fixes `log_s` for 1-parameter equations
   (a mazur fit's `opt$par` has no free `log_s`, and its `df`/AIC are unchanged from the pre-`log_s`
   build — a regression guard); (b) a 2-parameter fit estimates a free `log_s`; (c) **known-`(k, s,
   φ/σ_e)` recovery** on simulated GM and Rachlin data × both families (explicit relative-error asserts,
   geometric for `k`, looser tolerance for `s`); (d) **reduction-to-Mazur:** a GM (and Rachlin) fit on
   **Mazur-simulated data** estimates `s ≈ 1` within tolerance and recovers `k`, with `log_s` free — no
   user-fix mechanism is needed (the analytic `s = 1 ⇒ Mazur` reduction is already pinned in step 1).
   The test's direct `MakeADFun` (if any) also threads `map`.
4. **S3 + simulation** + `test-dd-tmb-methods.R`/`test-simulate-dd-mixed.R`: `s` appears in `tidy`/
   `summary`/`coef` only for 2-parameter fits; `predict`/`fitted` use the right mean; `simulate_dd_ip`
   recovers the input `s`.

## 6. Validation gates
- `R CMD check --as-cran` clean (still 1 env WARNING + 1 aspirational-URL NOTE); full suite 0 fail.
- Compile gate: C++ −nll == R density to 1e-8 for all four equations × both families.
- Recovery within tolerance on simulated GM/Rachlin data (both families); `s` recovered within ~0.15
  relative on adequately-informative designs.
- Reduction-to-Mazur at `s = 1` holds for GM and Rachlin.
- No regression: the existing mazur/exponential fits, S3 output, and emmeans are byte-compatible.

## 7. Risks
- **`pow(0, s)` in the AD tape** (Rachlin at `x = 0`) — **resolved by design** (§4.1): compute `pow` on a
  CondExp-guarded `x_safe` and select `mu = 1` for `x = 0`, so the tape never differentiates `pow` at 0.
  The compile gate includes a Rachlin `x = 0` row.
- **`s` identifiability** with few/clustered delays — kept population-level; the multi-start `log_s`
  perturbations + soft bounds mitigate; recovery tolerance is looser for `s` than `k`.
- **`map` threading** — `log_s` must be mapped on EVERY `MakeADFun` call (single + all multi-start **and
  the two tests that build `MakeADFun` directly**), or a 1-parameter fit would try to estimate an
  unidentified `log_s` (singular Hessian, wrong df/AIC). A test pins this.
- **Start sensitivity** — GM/Rachlin likelihoods can be flatter in `s`; the Mazur-based `k` start + `s = 1`
  start + multi-start should suffice, but watch exponential-family-style fragility (port the data-driven
  multi-start behavior).
- **Reporting placement of `s`** — a new `shape` row must not break broom-contract column counts in
  `tidy()` (exactly 8 columns) or downstream consumers; covered by a test.

## 8. Out of scope (later Phase-C slices)
Choice family (Family 2, Young 2018); subject-random `s` and/or `φ` (→ real `pdDiag`/`pdSymm`);
Ebert-Prelec (reduces to *exponential*, anchors to a different baseline); parametric draws;
per-group/covariate `s` design.
