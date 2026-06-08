# Design Spec — Choice family (Family 2): trial-level SS-vs-LL choice models

**Date:** 2026-06-08 · **Repo:** beezdiscounting · **Branch:** `feat/choice-family` (off `develop`)
**Phase:** Phase C.2 — the second Phase-C modeling slice (after the 2-parameter tier). Handoff
`docs/handoff-2026-06-07-phase-c-2param-tier.md` §2.1; SLT spec §11 "choice family".
**Scope:** add a trial-level binary-choice (smaller-sooner vs larger-later) family — a binomial GLMM via
TMB — with **two flavors in one template**: a **structural** model (subjective-value comparison via the
discount function; estimates `k` directly; shares the IP family's `k`/emmeans contract; enables the
IP↔choice tie-out) and a **descriptive** model (Young-2018 separate magnitude/delay sensitivities with
per-subject random slopes; no assumed discount function).

This spec was finalized after reading **Wileyto et al. (2004)** and **Young (2018)** directly; the model
choices below are deliberate syntheses of those validated forms with this package's architecture, not
defaults of convenience. See §1 for the reasoning.

## 1. Methodological anchors & the parameterization rationale
- **Wileyto, Audrain-McGovern, Epstein & Lerman (2004)**, *Behavior Research Methods* 36(1):41–51 —
  logistic regression recovers the hyperbolic `k`. Their Eq 3: `logit P(LL) = β1·(1 − 1/R) + β2·t`,
  `R = SS/LL`, `t = delay`, **no intercept**, `k = β2/β1` (Eq 4). Key lessons we use:
  - `k` as a **coefficient ratio** has an unstable interval — Wileyto found the **delta-method CI on
    `β2/β1` under-covered and fell back to a jackknife** (their Fig 4A vs 4B). → argues for estimating
    `k` *directly* with its own SE.
  - Their magnitude predictor is a **ratio** (`1 − 1/R = 1 − LL/SS`) → **scale-invariant** across reward
    magnitudes (a virtue we adopt).
  - `β1` is both the magnitude weight and the choice-curve steepness ("vagueness") — the two are
    **conflated**; there is no separate noise/sensitivity parameter.
  - No intercept, so `P = 0.5` sits exactly at the hyperbolic indifference — which is what makes `k`
    interpretable as "where subjective values balance."
- **Young, M. E. (2018)**, *JEAB* (doi:10.1002/jeab.316) — multilevel logistic regression treating
  **relative reward magnitude and relative delay as separate, per-subject-varying contributors**. His
  descriptive model (Eq 3, confirmed in Table 1): `logit P(LL) = β_mag·log(LL/SS) + β_delay·log(delay+1)`,
  **no intercept**, `β_mag > 0`, `β_delay < 0`, and **both weights vary randomly across subjects (random
  slopes)** — the central contribution. (For non-zero SS delay he uses `log((LL_delay+1)/(SS_delay+1))`;
  SS-delay is out of scope here.)
- **Our synthesis (the locked decisions):**
  - **Structural:** keep `k` an **explicit fitted parameter** (`k = exp(Xβ_k + σ_u u)`) — direct SE from
    `sdreport` (avoids Wileyto's ratio-CI problem), supports multiple discount functions, and inherits
    the IP family's factor/covariate/RE/emmeans machinery + the IP↔choice tie-out. Adopt Wileyto's
    **scale-invariance** via a **relative** value comparison `((ll/ss)·D(k,delay) − 1) = (V_LL−V_SS)/V_SS`
    (magnitude enters only as the `ll/ss` ratio). Keep a **separate sensitivity** `γ` (inverse
    temperature) — the construct Wileyto's `β1` conflates. The choice-**bias intercept `β0` is optional,
    default OFF** (off ⇒ `P=0.5` at indifference, `k` = the classical indifference rate, clean tie-out;
    on ⇒ models side bias at the cost of decoupling `k` from `P=0.5`).
  - **Descriptive:** adopt Young's exact predictors (`log(LL/SS)`, `log(delay+1)`, no intercept by
    default) and his **per-subject random slopes** (correlated 2-RE). This pulls the correlated-RE
    (Cholesky/LKJ) machinery into this slice for the *logistic* case.
- **Architecture reuse:** IP-family TMB machinery (`src/MixedDiscounting.h`, `R/dd-tmb.R`,
  `R/dd-tmb-methods.R`, `R/dd-comparisons.R`, `R/dd-param-space.R`): `.dd_discount_mu`/the `eqn_type`
  discount function, `build_fixed_rhs`/design + rank guard, the non-centered RE pattern, multi-start
  plumbing, param-space transforms, the `beta_k`-block emmeans contract. The 2-RE Cholesky pattern
  mirrors `beezdemand/src/HurdleDemand2RE.h` + `R/random-effects-utils.R`.

## 2. Scope
**In scope (this slice):**
- `fit_dd_choice()` — binomial GLMM, two `mode`s (structural | descriptive).
- **Structural**: equations **mazur + exponential** only (1-parameter discount functions; GM/Rachlin
  reuse the identical `.dd_discount_mu` and are a trivial follow-on — deferred to bound the slice). RE =
  **random intercept on `log k`**; between-subject factors/covariates on `log k` via `X`. Scale-invariant
  relative value comparison; optional `β0` (default off); separate `γ`.
- **Descriptive**: configurable fixed-effect formula on the choice logit (default `log(ll/ss) +
  log(delay+1)`, no intercept) with **correlated per-subject random slopes** on those predictors
  (2-RE Cholesky/LKJ) — matching Young 2018.
- New `src/ChoiceDiscounting.h` (one template, `mode` switch), `R/dd-choice.R`, `R/dd-choice-methods.R`,
  `R/simulate-dd-choice.R`, `mcq27_to_choice()` adapter, the 2-RE Cholesky helpers (logit case), tests,
  roxygen.
- The **IP↔choice tie-out** validation.

**Out of scope (later slices):** subject-random `γ`/`s`/`φ` in the *structural/discount* model (the
discount-parameter `pdSymm` case — distinct from the descriptive logistic random slopes done here);
GM/Rachlin in the structural value function; front-end (non-zero) SS delays; parametric-draw CIs;
random *intercept* + slopes (3-RE) variants of the descriptive model beyond the 2-slope default;
probability/effort-discounting choice.

## 3. Data contract
**Generic per-trial frame** (canonical input):
- `id` — subject id.
- `ss_amount` — smaller-sooner reward (> 0), available now (immediate).
- `ll_amount` — larger-later reward (> 0; expected `> ss_amount`, warned-not-aborted otherwise).
- `delay` — LL delay (> 0; SS delay fixed at 0 this slice).
- `choice` — `∈ {0,1}`, **1 = chose LL**.

**Validation** (`.dd_validate_choice`, mirroring `.dd_validate_ip`): coerce types; `choice` strictly
binary `{0,1}` (error, naming offending rows); `ss_amount`/`ll_amount` finite & > 0; `delay` finite & > 0;
one complete-case pass over all modeling columns (row-coherence); 0-indexed `subject_id` aligned to
**sorted** `subject_levels` (index by name, never positionally — the subject-id alignment trap). All
coercion/clamping is loud (warns + documents).

**MCQ27 adapter** `mcq27_to_choice(responses)`: joins trial-level MCQ27 responses
(`subjectid`/`questionid`/`response`) to the built-in 27-item lookup (each question's `ss_amount`,
`ll_amount`, `delay`) → the generic frame (`id=subjectid`, `choice=response`). Validates question ids
against the lookup. (Confirm the exact `mcq27`/`score_mcq27` item-table column names during planning —
the canonical Kirby 27-item parameter set.)

## 4. The model (one binomial template, `mode` switch)
Both modes: `choice_i ~ Bernoulli(p_i)`, `logit p_i = η_i`, non-centered subject random effects.

### 4.1 Structural mode (`mode = 0`)
- `D(k_i, delay_i)` = the `eqn_type` discount function (mazur `1/(1+k·x)`; exponential `exp(-k·x)`),
  reusing `.dd_discount_mu`. `V_SS_i = ss_amount_i` (immediate, `D(k,0)=1`); `V_LL_i = ll_amount_i ·
  D(k_i, delay_i)`.
- `k_i = exp( X_i·β_k + σ_u·u(subj_i) )` — the **same `log k = Xβ + σ_u u` predictor as the IP family**
  (RE = random intercept on `log k`; factors/covariates via `X`; the `k`/emmeans contract is shared).
- **Scale-invariant relative comparison:** `η_i = β0·[has_intercept] + γ·( (ll_amount_i/ss_amount_i) ·
  D(k_i,delay_i) − 1 )`, equivalently `γ·(V_LL−V_SS)/V_SS`. `γ = exp(log_gamma)` (population sensitivity,
  > 0); `β0` real (population choice-bias) **optional, default off** (`intercept = FALSE`).
- Free parameters: `β_k` (len `ncol(X)`), `log_sigma_u`, `log_gamma`, `beta0` (only if `intercept`),
  random `u` (scalar per subject).
- With `β0 = 0`, `η_i = 0` ⇔ `V_LL = V_SS` ⇔ `D(k,delay) = ss/ll` — the hyperbolic/exponential
  indifference; `k` is the classical indifference rate (this is the tie-out condition).

### 4.2 Descriptive mode (`mode = 1`)
- `η_i = Z_i·θ + Z^{re}_i·b_{subj_i}`, where `Z` is the fixed-effect design from a **configurable
  formula** (built via `build_fixed_rhs`/`model.matrix` + rank guard); `Z^{re}` is the random-effects
  design (default = the same magnitude/delay columns) and `b_subj ~ N(0, Σ)` are **correlated per-subject
  random slopes** (`Σ` a `q×q` covariance via a Cholesky factor — log-SDs + an LKJ/unconstrained
  correlation), mirroring `beezdemand`'s 2-RE pattern for the logit case.
- **Default formula:** `~ 0 + log(ll_amount / ss_amount) + log(delay + 1)` (no intercept; Young's exact
  predictors) → `η = θ_mag·log(ll/ss) + θ_delay·log(delay+1) + b_mag,subj·log(ll/ss) +
  b_delay,subj·log(delay+1)`. User-overridable via `predictors` (and an `re_formula`-style control for
  which terms get random slopes; default = both magnitude and delay).
- Free parameters: `θ` (len `ncol(Z)`), the Cholesky entries of `Σ` (log-SDs + correlations), random `b`
  (`q`-vector per subject). No `k`, `γ`, discount function, or scalar `u` (mapped out).

### 4.3 Shared
- Numerically stable binomial log-likelihood via the logit form (`dbinom_robust`/log-sum-exp); no
  `μ`-bound guard needed (logit link). The non-centered RE prior is `N(0, I)` on the standardized
  deviates; structural uses a scalar deviate, descriptive a `q`-vector transformed by the Cholesky `L`.

## 5. Backend `src/ChoiceDiscounting.h`
- `DATA_INTEGER(mode)` (0 structural, 1 descriptive); `DATA_INTEGER(eqn_type)` (structural: 0 mazur, 1
  exponential); `DATA_INTEGER(has_intercept)` (structural β0 on/off); `DATA_VECTOR(choice/ss_amount/
  ll_amount/delay)`, `DATA_IVECTOR(subject_id)`, `DATA_MATRIX(X)` (structural log-k design),
  `DATA_MATRIX(Z)` (descriptive fixed design), `DATA_MATRIX(Zre)` (descriptive RE design), `n_obs`,
  `n_subjects`, `n_re` (descriptive RE dim `q`).
- Parameters: `beta_k`, `log_sigma_u`, `log_gamma`, `beta0`, `theta`, `log_sd_re` (vector, len `q`),
  `cor_re` (unconstrained correlation params), `u` (scalar RE matrix, structural), `b` (`q`-col RE
  matrix, descriptive). Mode-irrelevant parameters are TMB **`map`-fixed** from R (structural maps
  `theta`/`log_sd_re`/`cor_re`/`b`; descriptive maps `beta_k`/`log_gamma`/`beta0`/`u`). Thread the `map`
  through **EVERY** `MakeADFun` call (single + multi-start + any direct test build) — the seam that bit
  the 2-parameter tier; one helper `.dd_choice_build_map(mode, has_intercept)` returns it; an
  extract-time assertion pins free-parameters ↔ mode.
- Structural η: `β0·has_intercept + γ·((ll/ss)·D(k,delay) − 1)`, `k = exp(X·β_k + σ_u·u(subj))`.
- Descriptive η: `Z·θ + (Zre ∘ (L·b(subj)))` with `L` the Cholesky factor of `Σ`.
- `nll -= dbinom_robust(choice_i, logit=η_i)` + `N(0,1)` prior on the standardized deviates.
- `ADREPORT` natural-scale `gamma`, `sigma_u`, and `Σ`/SDs+correlations; structural per-subject `k`
  reconstructed in R as in the IP family.
- Registered via `src/beezdiscounting_TMBExports.cpp` (`if (model == "ChoiceDiscounting") …`).

## 6. R API `R/dd-choice.R`
```
fit_dd_choice(data,
              mode = c("structural", "descriptive"),
              id_var = "id", ss_var = "ss_amount", ll_var = "ll_amount",
              delay_var = "delay", choice_var = "choice",
              equation = c("mazur", "exponential"),     # structural
              intercept = FALSE,                          # structural β0 (default off)
              factors = NULL, factor_interaction = FALSE, # structural: on log k
              continuous_covariates = NULL,
              predictors = NULL,                          # descriptive RHS;
                                                          #   default ~ 0 + log(ll/ss) + log(delay+1)
              random_slopes = TRUE,                       # descriptive: per-subject slopes (2-RE)
              start_values = NULL, tmb_control = ..., multi_start = TRUE, verbose = 1, ...)
```
- validate → prepare (0-indexed subjects, row-coherence) → build `X` (structural) or `Z`/`Zre`
  (descriptive) → assemble TMB data + the mode/intercept-appropriate `map` → multi-start optimize →
  `sdreport` → assemble S3 `beezdiscounting_choice` (parallels `beezdiscounting_tmb`: `call`, `opt`,
  `model$coefficients/se`, `sdr`, `param_info` incl. `mode`/`equation`/`intercept`/factor metadata/
  RE spec, `formula_details`, `subject_pars` (structural: id/u_i/k), `loglik/AIC/BIC`, `converged`,
  `se_available`, `data`, `coercion_info`).
- `mcq27_to_choice()` adapter (§3).

## 7. S3 surface `R/dd-choice-methods.R`
`tidy / glance / summary / confint / coef / fixef / ranef / VarCorr / predict / fitted / residuals /
augment / logLik / AIC / BIC / nobs / print` for `beezdiscounting_choice`, mirroring the IP contracts and
reusing `.dd_transform_coef_table`/param-space:
- **predict**: `type = "prob"` → fitted `P(LL)` (subject and/or population level, nlme-style columns);
  `type = "parameters"` → subject `k` table (structural).
- **Structural** fits expose `k` on the `log k` scale and route through `get_dd_param_emms`/
  `get_dd_comparisons` for the **shared `k`/emmeans contract** (verify the `beta_k`-block covariance
  selection still isolates with the extra `log_gamma`/`beta0` present — keep `beta_k` first in the param
  order, the property proven for `log_s`; add an isolation guard test). `γ` (and `β0` if on) surface as
  their own rows (`γ = exp(log_gamma)` via param-space).
- **Descriptive** fixed sensitivities (`θ`) via `tidy`/`summary`; the random-slope covariance via
  `VarCorr`/`ranef` (per-subject magnitude/delay sensitivities — Young's deliverable). No dedicated
  emmeans this slice.
- `glance` reports `mode`, `equation`/`intercept` (structural), backend, nobs, n_subjects, n_re,
  converged, logLik/AIC/BIC.

## 8. Simulation `R/simulate-dd-choice.R`
`simulate_dd_choice(n_subjects, design | n_trials, mode, equation, log_k_pop, sigma_u, gamma, beta0,
theta, re_sd, re_cor, seed, …)`:
- **structural**: `k_i = exp(log_k_pop + σ_u u_i)`; `p = plogis(β0 + γ((ll/ss)·D(k,delay) − 1))`; draw
  `choice ~ Bernoulli(p)`.
- **descriptive**: draw correlated per-subject slopes `b_i ~ N(0, Σ)`; `p = plogis(Z θ + Zre b_i)`; draw
  choices.
- Returns the generic per-trial frame. Known-truth generator for recovery + the tie-out.

## 9. The IP↔choice tie-out (key correctness test)
Simulate **IP** data and **structural choice** data from the **same** `(log_k_pop, σ_u, equation)` with
`β0 = 0` (default). Fit `fit_dd_tmb` (IP) and `fit_dd_choice(mode="structural")`; assert recovered
population `k = exp(β_k[1])` agrees within tolerance (geometric; looser than within-family recovery).
Pins that the structural choice model measures the *same* discount construct as the IP model. (The
scale-invariant relative comparison makes the indifference `V_LL=V_SS` coincide with the IP indifference
exactly when `β0=0`.)

## 10. Build order (TDD — test named first)
1. **R reference + value/probability helpers** (pure R, no compile): `V`/relative-comparison via
   `.dd_discount_mu`; structural `plogis(β0 + γ((ll/ss)D − 1))`; descriptive design from the default
   formula + a Cholesky→Σ helper. Assert: scale-invariance (same `ll/ss`, `delay` ⇒ same `p`
   irrespective of absolute magnitude); structural `β0=0` ⇒ `p=0.5` at `V_LL=V_SS`; reduction sanity.
2. **C++ template + compile gate**: `src/ChoiceDiscounting.h`; gate asserts C++ −nll == R binomial nll to
   1e-8 for **both modes** (structural mazur/exp, β0 on & off; descriptive default design with a fixed
   Cholesky), `map` threaded per mode. Include `choice ∈ {0,1}` rows and a near-saturated `p`.
3. **2-RE Cholesky machinery (descriptive)** + tests: log-SDs + correlation → `Σ` round-trip; the
   non-centered `L·b` transform; a known-`Σ` recovery on simulated descriptive data.
4. **Fit pipeline + map**: `fit_dd_choice` validate→prepare→design→map→optimize→extract; per-mode `map`
   correctness (structural: no free `θ`/RE-cov; descriptive: no free `β_k`/`γ`/`β0`/`u`), df/AIC,
   multi-start guard, the `intercept` on/off path.
5. **Recovery + tie-out**: structural `(k, γ[, β0])` recovery; descriptive `(θ, Σ)` recovery; the
   IP↔choice tie-out (§9). Explicit relative-error asserts; geometric for `k`; **centered seeds**
   (per the 2-param seed-robustness lesson — a noisy cell gets a seed sweep + a centered seed, never a
   cherry-picked pass).
6. **S3 + simulation + MCQ27 adapter**: method contracts; `predict` prob; structural emmeans isolation
   guard; `VarCorr`/`ranef` for the descriptive slopes; `mcq27_to_choice` round-trip; `simulate_dd_choice`
   recovers inputs in expectation.

## 11. Validation gates
- `R CMD check --as-cran` clean except the known env-clang WARNING + aspirational-Pages-URL NOTE (no new
  findings).
- Compile gate: C++ −nll == R binomial nll to 1e-8 (both modes).
- Recovery within tolerance (structural `k`/`γ`[/`β0`]; descriptive `θ`/`Σ`); IP↔choice tie-out holds.
- No regression to the IP family / `get_dd_param_emms` / the 2-parameter tier (full suite green).
- **Codex CLI review of the plan before coding AND a Codex CLI capstone review of the implemented feature
  before merge** (both via the CLI, read-only — the standing convention).

## 12. Risks
- **Descriptive 2-RE identifiability / convergence** — correlated random slopes on two predictors need
  adequate per-subject trials and stimulus-space coverage; multicollinearity of magnitude & delay (Young
  Fig 1) can destabilize. Mitigate: non-centered Cholesky parameterization, sensible SD/correlation
  starts + soft bounds, multi-start; document data requirements. This is the newest machinery — port
  carefully from `beezdemand/src/HurdleDemand2RE.h` and verify the `Σ` round-trip in isolation (step 3).
- **`map` threading across modes** (structural vs descriptive fix different parameter blocks) — one
  helper + an extract-time free-param↔mode assertion; thread through EVERY `MakeADFun`.
- **emmeans isolation** with extra structural parameters (`log_gamma`, optional `beta0`) in `opt$par`/
  `cov.fixed` — keep `beta_k` first so the `beta_k`-block selection is unshifted (proven for `log_s`);
  isolation guard test.
- **Binomial separation** on near-deterministic choosers (all-SS/all-LL) — the RE + soft bound/sensible
  start on `log_gamma` (structural) and shrinkage from the random slopes (descriptive) mitigate; mirror
  the IP φ-floor pattern if a degenerate sink appears.
- **Structural magnitude scaling** — resolved by the scale-invariant relative comparison (§4.1); a
  pure-R test asserts magnitude scale-invariance.
- **Descriptive transform domains** — `log(delay+1)` needs `delay ≥ 0` (validated; `+1` handles 0);
  `log(ll/ss)` needs both amounts > 0 (validated). The configurable formula lets users choose other
  transforms (note the default = Young's exact form).

## 13. Out of scope (later slices)
Subject-random `γ`/`s`/`φ` in the structural/discount model (discount-parameter `pdSymm`); GM/Rachlin in
the structural value function; front-end (non-zero) SS delays; parametric draws; a descriptive-mode
emmeans/contrast layer; 3-RE (intercept + slopes) descriptive variants; probability/effort-discounting.
