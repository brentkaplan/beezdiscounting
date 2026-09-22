# beezdiscounting 0.4.0

This is a large release (the first since 0.3.2, January 2025): mixed-effects
and Bayesian discounting tiers, trial-level choice models, 21-item MCQ and PDQ
scoring, Monte Carlo power analysis, and ten vignettes (all new since 0.3.2,
which shipped none).

### Bug fixes that can change results

* `check_unsystematic()` now applies both Johnson & Bickel (2008) thresholds
  on the `y` scale (`c1 * ll`, `c2 * ll`). Previously `ll` cancelled out, so
  amount-scale data (e.g. `y` in dollars with `ll = 100`) were judged against
  the bare proportions `c1`/`c2` and received wrong verdicts. Criterion 2 is
  now strict: a decline of exactly `c2 * ll` passes. `ll`, `c1` and `c2` are
  validated. Verdicts on proportion data with `ll = 1` change only for
  changes of exactly `c1` or `c2`, which floating-point error previously
  judged either way (e.g. `y = c(0.6, 0.8, 0.1)` used to fail C1).
* `check_unsystematic()` and `calc_aucs()` now error on missing `y`/`x`
  values and on duplicate delays within a subject, which previously gave a
  silent `NA`/`TRUE` verdict or a result that depended on row order.
  `calc_aucs()` documents that `auc_log10` uses `log10(x + 1)` and therefore
  depends on the delay unit.
* INN imputation (`score_mcq()`, `score_mcq27()`, `score_pdq()` with
  `impute_method = "inn"`) now warns when items stay missing (neighbours
  disagree with `random = FALSE`, or a whole rank group is missing) instead
  of letting the resulting `NA` scores pass silently.
* `score_dd()`, `ans_dd()` and `calc_dd()` read numeric Qualtrics exports per
  item. The bundled 5.5-trial template codes `Attend-LL` in reverse order
  (1 = "in 25 years", 2 = "now"), so numeric exports previously flagged
  respondents who chose "now" and passed those who chose the 25-year option,
  and applied the attention-item k to the wrong respondents. Text exports were
  unaffected.
* `score_pd()`, `ans_pd()` and `calc_pd()` have the same fix: the 5.5-trial
  probability template codes `Attend-LL` in reverse order (1 = "with a 1%
  chance", 2 = "for sure"), so numeric exports previously flagged respondents
  who chose the certain amount and passed those who chose the 1% chance.
  Text exports were unaffected.
* The 5.5-trial delay and probability scorers (`score_dd()`, `ans_dd()`,
  `calc_dd()`, `score_pd()`, `ans_pd()`, `calc_pd()`) now recognise only the
  template's numeric codes and option texts. Previously any other value
  (a typo, `"-99"`, a relabelled export, a blank cell) was silently scored as
  the delayed/uncertain choice with a real k or h. Blank cells are now
  treated as unanswered and dropped, and other unrecognised values raise an
  error naming the item and value.
* `fit_dd_tmb()` and `fit_dd_choice()` (both modes) no longer let a
  non-converged start displace a converged one on negative log-likelihood
  alone: the multi-start keeps the lowest-NLL converged start that passes the
  log-k sanity guard, falling back to other starts only when no start both
  converged and passed the guard (recorded in the new `multi_start_info`
  element). A non-converged
  fit now raises a classed `beezdiscounting_convergence_warning` at fit time
  regardless of `verbose`, and `tidy()`, `confint()`, `summary()`,
  `get_dd_param_emms()` and `get_dd_comparisons()` repeat it. Previously such
  fits (e.g. `tmb_control = list(iter_max = 3)`) returned finite p-values and
  intervals with no warning. Non-positive-definite Hessian and failed
  `sdreport()` warnings are likewise no longer silenced by `verbose = 0`.
* `simulate_dd_ip(family = "gaussian")` no longer clamps draws to `[0, 1]`,
  and `fit_dd_tmb(family = "gaussian")` / `fit_dd_brms(family = "gaussian")`
  keep responses as supplied: no clamping to `[0, 1]` and no automatic
  percent detection (data that look percent-scaled get a warning; an explicit
  `response_scale = "percent"` or `"amount"` still converts). The Gaussian likelihood is unbounded, so the clamp made
  the simulated data, and the data the model saw, differ from the model
  being fitted: near delay 0 about half the draws piled up at 1.
  Gaussian-family simulations and `power_discounting(family = "gaussian")`
  results change; the SLT-beta family is unaffected.

### Monte Carlo power analysis

* `power_discounting()` estimates statistical power for detecting a
  between-subject difference in discount rate (`delta_k` on log k) by
  simulating with `simulate_dd_ip()` and refitting each replicate with
  `fit_dd_tmb()`. Reports the power estimate with a Wilson Monte Carlo
  confidence interval, p-value and CI-exclusion hit rates, and convergence
  diagnostics; non-usable fits are excluded from the denominator and
  surfaced, never counted as misses. The Wald test uses a t reference with
  the design's two-sample df (`n - 2`): an empirical small-sample
  calibration, not a model-derived df, chosen because it keeps the Type I
  error near nominal in the package's calibration tests where the z-test
  runs high.
* `find_n_discounting()` searches for the smallest total N reaching a
  target power via bisection, adding replicates adaptively where the Monte
  Carlo verdict is ambiguous and re-confirming the selected N before
  reporting.
* Type I error calibration, convergence handling, a closed-form benchmark
  against `pwr::pwr.t.test()`, monotonicity, and seed reproducibility are
  verified in the test suite; see `vignette("power-analysis")` for scope
  and validity notes. Mirrors `beezdemand::power_demand()`; the within- vs
  between-subject asymmetry is intentional.

### Probability Discounting Questionnaire (PDQ)

* New `score_pdq()` scores the 30-item PDQ (Madden, Petry, & Johnson,
  2009): per-block h under the hyperbolic odds-against model, consistency,
  and risky choice ratios, with the same strict input validation and
  imputation options as `score_mcq()`. The implementation reproduces the
  Gray et al. (2016) scoring-syntax lookup tables exactly for all 3 x 1024
  response patterns (verified in the test suite). An `overall_h` from a
  pooled 30-item ladder is also reported as a documented beezdiscounting
  extension (the published scoring defines no overall ladder; `mean_h` is
  Gray et al.'s recommended composite).
* New `prop_sc()` (guaranteed-choice proportions by h rank) and
  `pdq_to_choice()` (trial-level choice frame including the odds against
  winning, `theta`), plus a bundled `pdq` example dataset and a
  "Scoring the Probability Discounting Questionnaire" vignette.
* `get_lookup_table()` gains an `instrument` argument
  (`"mcq27"`, `"mcq21"`, `"pdq"`); `items` remains as a back-compatible
  alias. Internally the MCQ registry is now instrument-keyed and the
  ladder-scoring core is shared across instruments; 27- and 21-item MCQ
  results are unchanged (pinned by golden-fixture regression tests).

### 21-item MCQ support

* New `score_mcq()` scores both the 27-item (Kirby, Petry, & Bickel, 1999)
  and the original 21-item (Kirby & Marakovic, 1996) MCQ with the same
  consistency-maximization algorithm; `score_mcq27()` is now a wrapper and
  its output is unchanged. The 21-item item table (amounts, delays, k at
  indifference, ranks) was transcribed from Kirby & Marakovic (1996,
  Table 1) and cross-validated against the Kaplan et al. (2014) Excel
  Automated Scorer, including its ladder-edge conventions.
* `prop_ss()`, INN missing-data imputation (now defined over k-rank
  neighbor groups, identical results for the 27-item version), the new
  `mcq_to_choice()` (generalizing `mcq27_to_choice()`), and
  `get_lookup_table()` all support `items = 21`.
* New `mcq21` example dataset.
* Stricter input validation in the scorer: duplicate, unknown, missing, or
  fractional/non-whole question ids and responses outside 0/1/NA now error
  instead of silently producing invalid scores. `mcq_to_choice()` still
  accepts ragged input.
* Character, factor, and logical responses to `score_mcq()`/`score_mcq27()`
  are now normalized and score identically to numeric 0/1 (previously they
  either passed validation and then failed obscurely during scoring, or
  errored).
* New `plot()` method for 21-item `score_mcq()` output
  (`plot.score_mcq_output()`), matching the existing 27-item plot method.
* Internal (unexported) `inn()` gained a `reg` parameter --
  `inn(dat, reg, random, verbose)` -- to support both MCQ versions. This is
  a breaking signature change for any code calling
  `beezdiscounting:::inn()` directly.


### New vignettes

* `vignette("mcq27-scoring")`: scoring the 27-item Monetary Choice Questionnaire
  (`score_mcq27()`, `get_lookup_table()`, `prop_ss()`, `summarize_mcq()`,
  `mcq27_to_choice()`, and the wide/long converters).
* `vignette("delay-discounting-basics")`: a getting-started walk-through of
  indifference-point screening (`check_unsystematic()`), curve fitting
  (`fit_dd()`/`results_dd()`), discount rate `k`, and AUC (`calc_aucs()`).
* `vignette("fivetrial-task")`: scoring the 5.5-trial delay and probability
  discounting tasks from the Qualtrics minute-discounting template.
* `vignette("tmb-mixed-effects")`: the `fit_dd_tmb()` mixed-effects workflow and
  its S3 methods, prediction, diagnostics, and group comparisons.

### Bug fixes

* `tidy()`, `confint()`, and `summary()` on a `fit_dd_tmb()` or structural
  `fit_dd_choice()` fit with `factors` reported the intercept's standard
  error for every `beta_k` coefficient (the log-k intercept and each
  condition contrast): the internal SE lookup collapsed the duplicated
  `beta_k` parameter names onto the first element. Condition-contrast Wald
  statistics, p-values, and confidence intervals were therefore
  anticonservative (the `power_discounting()` Type I calibration battery
  measured a false-positive rate of 0.168 at nominal .05 before the fix).
  SEs are now aligned positionally with the coefficient vector in both
  accessors, and ambiguous legacy objects with duplicated names return `NA`
  SEs instead of silently misaligned values.

* `check_unsystematic()` and `calc_aucs()` now compute results per `id`. Previously
  they computed a single result over the whole data frame and recycled it across
  `unique(id)`, so multi-subject input returned the same verdict / AUC for every
  subject. Each now returns one correct row per subject (single-subject output is
  unchanged). `check_unsystematic()` orders points by `x` when that column is present,
  and `calc_aucs()` orders by `x` within each subject. Rows with a missing `id` are
  dropped so they cannot contaminate other subjects' results.
* `prop_ss()` now pools correctly across respondents. It previously dropped all but
  the first occurrence of each `questionid` (via `match()`) and divided by a fixed
  denominator of 3, so multi-respondent input was reduced to one respondent and could
  return values above 1. It now averages over all retained rows per k-rank and returns
  `NA` for a rank with no observed responses.

### Subject-random `s` (GM/Rachlin curvature)

* `fit_dd_tmb(..., random_effects = k + s ~ 1)` now fits a per-subject random
  intercept on the Green-Myerson / Rachlin curvature `s`, jointly with `log k`,
  with `covariance_structure = "pdSymm"` (correlated) or `"pdDiag"`
  (independent), for both `family = "sltb"` and `"gaussian"`. Per-subject `s_i`
  is soft-clamped toward `(0.05, 20)` (a C-infinity softplus map; see the de-hang
  note below). `VarCorr()`/`ranef()`/`summary()` surface the
  `(k, s)` covariance and per-subject `s`. `predict()` at
  `level = "subject"` uses each subject's estimated `s_i`; response-SD
  standardization uses the population precision `phi`. `simulate_dd_ip()` gains
  `sigma_s`/`rho_ks` to generate `(log k, log s)` bivariate-normal data for
  recovery testing.
* The per-subject `s` random effect now uses a smooth (C-infinity) softplus soft
  clamp instead of a hard clamp, so a clamp-binding subject (e.g. a no-discounting
  or step subject) converges instead of hanging the Laplace inner solve.

### Bayesian (brms) modeling tier

* New `fit_dd_brms()`: Bayesian mixed-effects discounting via brms/Stan for
  all four TMB equations (`"mazur"`, `"exponential"`, `"green-myerson"`,
  `"rachlin"`). `family = "beta"` (identity link with a differentiable
  squish) is an ordinary beta likelihood, not the TMB SLT-beta, and its
  estimates are not expected to match `fit_dd_tmb(family = "sltb")` when
  responses sit at or near 0 or 1; `family = "gaussian"` matches the TMB
  gaussian likelihood wherever the TMB mu clamp does not bind (everywhere
  except extreme decay underflow). Exact 0/1 responses must be handled
  explicitly: the default `boundary = "error"` refuses to fit and reports
  the count; `"squeeze"` (Smithson-Verkuilen, rescales every response) and
  `"zoib"` (zero-one-inflated beta; changes the estimand) are opt-ins.
  `summary()` reports the exact- and near-boundary fractions.
  `init = "tmb"` (both brms fitters) seeds the chains from the TMB pre-fit
  only when that fit converged with finite estimates, and otherwise falls
  back to prior-center inits with a warning.
* New `fit_dd_choice_brms()`: the structural trial-level choice model under
  `bernoulli("logit")`, matching `fit_dd_choice(mode = "structural")`,
  including between-subject designs on `log k` via
  `factors`/`factor_interaction`/`continuous_covariates` (rank-deficient
  designs are rejected before sampling; `gamma` and `b0` stay
  population-level).
* New `default_dd_priors()` / `default_dd_choice_priors()`: inspectable
  defaults with delay-unit-aware `logk` anchoring (centers
  `k * median(delay) = 1`).
* `fit_dd_brms(random_effects = k + phi ~ 1, family = "beta")` adds a
  per-subject precision random effect, correlated with the `log k` intercept
  (`covariance_structure = "pdSymm"`, the default) or independent
  (`"pdDiag"`). `phi` becomes a predicted distributional parameter on the log
  link; `variance_components` and `subject_pars` gain the precision-RE SD,
  the `(k, phi)` correlation, and per-subject precision, mirroring
  `fit_dd_tmb(random_effects = k + phi ~ 1)`. The TMB per-subject precision
  floor (0.1) is not replicated by the brms log-link RE.
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
* New vignette `vignette("bayesian-discounting")` (precomputed output) walks
  through `fit_dd_brms()` / `fit_dd_choice_brms()` and their S3 surface.
* New vignette "Comparing discounting rates between groups"
  (`vignette("dd-group-comparisons")`): factor designs on log k, estimated
  marginal means, and contrasts across both backends -- TMB (Wald + holm)
  and brms (posterior draws + `post.prob`) -- for indifference-point and
  trial-level choice models alike.

### Modeling tiers and choice models

* New `fit_dd_choice(mode = "structural")` fits trial-level smaller-sooner vs
  larger-later choice as a binomial GLMM, estimating the discount rate `k`
  directly (scale-invariant value comparison, optional choice-bias intercept).
  It shares the `get_dd_param_emms()`/`get_dd_comparisons()` `k` contract with
  `fit_dd_tmb()` and is validated by an IP-vs-choice tie-out. `simulate_dd_choice()`
  generates structural choice data.

* `fit_dd_choice(mode = "descriptive")` fits the Young (2018) descriptive
  choice model: a logistic mixed model on the log amount ratio and log delay
  with subject random slopes, returning logit-scale sensitivities rather than
  a discount rate. Shares the S3 surface (`tidy()`, `summary()`, `predict()`,
  ...) with the structural mode; see `vignette("choice-discounting")`.

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

* New `plot_qq()` methods (re-exporting `beezdemand::plot_qq()`) for
  `fit_dd_tmb()` and `fit_dd_choice()` fits: a normal QQ plot of the
  estimated (shrunken, empirical-Bayes) subject random-effect deviates
  against a normal reference, the standard check on the Gaussian
  random-effects assumption. Bayesian (`fit_dd_brms()`) fits are
  intentionally excluded; use `brms::pp_check()` and MCMC diagnostics
  there instead.

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

### Documentation (TMB tier)

- New vignette `sltb-discounting`: why bounded error distributions matter for
  indifference points, the SLT-beta density, a boundary demonstration on the
  bundled example data, and the mixed-effects workflow.

### Notes

- The data validator coerces percent/amount response scales to `[0, 1]` and,
  for the SLT-beta and beta families, clamps mild out-of-range values,
  **warning loudly** and naming the number of values coerced or clamped
  (Gaussian fits keep responses as supplied).

### Bug fixes (scoring)

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
