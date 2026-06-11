# beezdiscounting (development version)

### Bayesian (brms) modeling tier

* New `fit_dd_brms()`: Bayesian mixed-effects discounting via brms/Stan for
  all four TMB equations (`"mazur"`, `"exponential"`, `"green-myerson"`,
  `"rachlin"`). `family = "beta"` (identity link with a differentiable
  squish) is the closest brms analog of the TMB SLT-beta;
  `family = "gaussian"` matches the TMB gaussian likelihood wherever
  the TMB mu clamp does not bind (everywhere except extreme decay
  underflow). Boundary responses are
  handled via Smithson-Verkuilen squeezing (default), zero-one-inflated
  beta (`boundary = "zoib"`; changes the estimand), or refusal.
* New `fit_dd_choice_brms()`: the structural trial-level choice model under
  `bernoulli("logit")`, matching `fit_dd_choice(mode = "structural")`,
  including between-subject designs on `log k` via
  `factors`/`factor_interaction`/`continuous_covariates` (rank-deficient
  designs are rejected before sampling; `gamma` and `b0` stay
  population-level).
* New `default_dd_priors()` / `default_dd_choice_priors()`: inspectable
  defaults with delay-unit-aware `logk` anchoring (centers
  `k * median(delay) = 1`).
* S3 methods mirror the TMB tier's contracts: the exact 8-column `tidy()`
  table (with `NA` `statistic`/`p.value`; estimates are posterior medians of
  report-space-transformed draws), `glance()` with `elpd_loo`/`looic` and
  MCMC diagnostics in place of `AIC`/`BIC`, `confint()` quantile credible
  intervals, `predict()`, `ranef()` with per-subject `k` posterior
  summaries for indifference-point fits, `print()`/`summary()`.
* `get_dd_param_emms()` and `get_dd_comparisons()` accept brms
  indifference-point and structural choice fits: draws-based marginal
  means and contrasts over the same reference grid as the TMB path, with
  quantile credible intervals and `post.prob` (posterior probability of
  direction) in place of adjusted p-values.
* `init = "tmb"` is available in both Bayesian fitters (a quiet TMB
  pre-fit supplies the chain starting values, with prior-center fallback).
* brms, posterior, and loo are Suggests-only.

# beezdiscounting 0.4.0

### New Features

* New `fit_dd_choice(mode = "structural")` fits trial-level smaller-sooner vs
  larger-later choice as a binomial GLMM, estimating the discount rate `k`
  directly (scale-invariant value comparison, optional choice-bias intercept).
  It shares the `get_dd_param_emms()`/`get_dd_comparisons()` `k` contract with
  `fit_dd_tmb()` and is validated by an IP-vs-choice tie-out. `simulate_dd_choice()`
  generates structural choice data. (Descriptive Young-2018 model: forthcoming.)

* New `mcq27_to_choice()` reshapes long-form 27-item Monetary Choice
  Questionnaire responses (`subjectid`/`questionid`/`response`) into the
  per-trial `id`/`ss_amount`/`ll_amount`/`delay`/`choice` frame consumed by
  `fit_dd_choice()`, using the canonical Kirby, Petry, & Bickel (1999) item
  design. `get_lookup_table()` now returns that complete design (adds
  `ss_amount`, `ll_amount`, and `delay`).

* `fit_dd_tmb()` and `simulate_dd_ip()` gain the two-parameter hyperboloid
  equations `"green-myerson"` (`mu = (1 + k*x)^(-s)`) and `"rachlin"`
  (`mu = 1 / (1 + k*x^s)`), with a single population nonlinearity exponent `s`
  (reported by `tidy()`/`summary()`/`confint()`). Both reduce to `"mazur"` at
  `s = 1`.

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

### Documentation

- New vignette `sltb-discounting`: why bounded error distributions matter for
  indifference points, the SLT-beta density, a boundary demonstration on the
  bundled example data, and the mixed-effects workflow.

### Notes

- The data validator coerces percent/amount response scales to `[0, 1]` and
  clamps mild out-of-range values, **warning loudly** and naming the number of
  values coerced or clamped.

### Bug Fixes

- **`score_dd()`** and **`ans_dd()`**: Fixed incorrect response classification for
  Qualtrics numeric recode exports where SS = `"1"` and LL = `"2"`. Previously,
  all numeric responses were classified as `"ll"`, producing incorrect `kval` and
  `ed50`. Both text exports (containing `"now"`) and numeric exports (`"1"` / `"2"`)
  are now handled correctly via an internal `normalize_dd_response()` helper.

- **`score_pd()`** and **`ans_pd()`**: Same fix for probability discounting.
  SC = `"1"`, LU = `"2"` in numeric exports are now correctly classified via
  an internal `normalize_pd_response()` helper.

# beezdiscounting 0.3.2

### New Features

- **`fit_dd()`**:
  - Introduced a new function to fit delay-discounting models using specified equations (`"mazur"`/`"hyperbolic"` or `"exponential"`) and methods (`"pooled"`, `"mean"`, or `"two stage"`).
  - Supports flexible data handling for aggregated and participant-specific modeling.
  - Returns an object of class `"fit_dd"` containing the fitted models, input data, and method details.

- **`plot_dd()`**:
  - Added a function to visualize fitted delay-discounting models.
  - Automatically adapts to different fitting methods, including aggregated and individual models.
  - Provides customizable axis labels, title, and optional log-transformed x-axis for improved visualization of delay scales.

- **`results_dd()`**:
  - New utility to extract model parameter estimates, confidence intervals, and fit statistics from a `"fit_dd"` object.
  - Supports both aggregated and participant-specific models.
  - Outputs a tidy tibble with columns for terms, estimates, standard errors,
    t-statistics, p-values, R2, three different AUC metrics, and confidence bounds.

- **`check_unsystematic()`**:
  - New utility function to check delay-discounting datasets for unsystematic
    data patterns according to Johnson & Bickel's (2008) two criteria.

- **`calc_aucs()`**:
  - New utility function to calculate three different area under the curve
    (AUC) metrics for delay-discounting data according to Borges et al. (2016).

### Improvements

- Confidence intervals are now computed using the `calc_conf_int()` function, ensuring accurate estimation based on model degrees of freedom.
- R2 values are calculated consistently using the `calc_r2()` function, providing reliable fit metrics for all models.

### Enhancements

- The package now supports robust delay-discounting workflows, from unsystematic
  identification (`check_unsystematic`), model fitting (`fit_dd`), to visualization (`plot_dd`), to result extraction
  (`results_dd`).
- Improved compatibility with delay-discounting datasets that require participant-level or aggregated modeling approaches.


# beezdiscounting 0.3.1

## Minor fix

* Correctly names output columns from `calc_pd()` and `score_pd()`. `ep50` changed to `etheta50` and corrected calculation of `ep50`.

# beezdiscounting 0.3.0

## New features

* Add functions for scoring 5.5 trial probability discounting task (from the Qualtrics template) including: `calc_pd()`
(and `score_pd()`, `timing_pd()`, and `ans_pd`).

## Minor fix

* Subsetting issue is fixed in `score_dd()` that would unintentionally drop all rows if both conditions were `FALSE`.

## Other changes

* Rename example data from `five.fivetrial` to `five.fivetrial_dd` for delay discounting.

* Add example data `five.fivetrial_pd` for probability discounting.

# beezdiscounting 0.2.0

## New features

* `score_mcq27()` properly supports arguments: `impute_method`, `random`, `return_data`, and `verbose`.
See documentation and the `README` for explanations.

* `generate_data_mcq()` can generate fake MCQ data, including `seed` and `prop_na` arguments for
reproducibility and specifying proportion of `NA`s.

* `long_to_wide*` and `wide_to_long*` are helper functions to reshape data from/to different formats.

## Minor fix

* When no imputation is specified and `NA`s exist in the data, `score_mcq27()` returns `NA`s for the scoring
instead of 1.

# beezdiscounting 0.1.0

* Initial release with basic scoring of 27-item Monetary Choice Questionnaire and 5.5 trial delay discounting task from the Qualtrics template.

* Added a `NEWS.md` file to track changes to the package.
