# Fit an indifference-point mixed-effects discounting model via TMB

Fits a discounting model (Mazur hyperbolic, exponential, or the
two-parameter Green-Myerson / Rachlin hyperboloids) with a random
intercept on `log k`, between-subject fixed effects, and either an
SLT-beta or Gaussian observation family, using Template Model Builder
for exact AD + Laplace approximation.

## Usage

``` r
fit_dd_tmb(
  data,
  y_var = "y",
  x_var = "x",
  id_var = "id",
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
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
  covariance_structure = c("pdSymm", "pdDiag"),
  ...
)
```

## Arguments

- data:

  Long data frame with subject id, delay, and indifference proportion
  columns.

- y_var, x_var, id_var:

  Column names (defaults `"y"`, `"x"`, `"id"`).

- equation:

  One of `"mazur"`, `"exponential"`, `"green-myerson"`, or `"rachlin"`.
  The two 2-parameter (hyperboloid) forms add a single population
  nonlinearity exponent `s` (estimated on the log scale) and reduce to
  `"mazur"` at `s = 1`.

- family:

  Observation family: `"sltb"` (default) or `"gaussian"`. For `"sltb"`,
  responses outside `[0, 1]` after scaling are clamped (with a warning);
  for `"gaussian"`, whose likelihood is unbounded, they are kept as
  observed and percent-scaled data are not detected automatically (set
  `response_scale = "percent"`).

- random_effects:

  RE formula: `k ~ 1` (single random intercept on `log k`),
  `k + phi ~ 1` (a joint 2-D random intercept on `(log k, log phi)`,
  SLT-beta only), or `k + s ~ 1` (a joint 2-D random intercept on
  `(log k, log s)`, Green-Myerson and Rachlin only).

- factors:

  Character vector of between-subject factor names.

- factor_interaction:

  Logical; include a pairwise factor interaction.

- continuous_covariates:

  Character vector of covariate names.

- ll:

  Optional larger-later reward for `amount`-scale coercion.

- response_scale:

  One of `"proportion"`, `"percent"`, `"amount"`.

- start_values:

  Optional named list overriding defaults.

- tmb_control:

  Optimizer control list.

- multi_start:

  Logical; if `TRUE` (default), run the 3-set guarded multi-start.
  Converged starts are preferred over non-converged ones, then the
  lowest negative log-likelihood wins; the Hessian is checked on the
  kept fit only.

- verbose:

  Integer verbosity (0 silent, 1 progress, 2 debug).

- covariance_structure:

  Covariance for a 2-D random effect (`k + phi ~ 1` or `k + s ~ 1`):
  `"pdSymm"` (default; correlated random intercepts) or `"pdDiag"`
  (independent, correlation fixed at 0). Ignored for `k ~ 1`.

- ...:

  Reserved.

## Value

An object of class `beezdiscounting_tmb` with components:

- call:

  The matched call.

- opt:

  Normalized optimizer result (`par`, `objective`, `convergence`,
  `message`).

- model:

  List of `coefficients`, `se`, and `variance_components`.

- sdr:

  TMB `sdreport` object (or `NULL` if SE computation failed).

- hessian_pd:

  Logical positive-definiteness of the Hessian.

- param_info:

  Model metadata (equation, family, dimensions, factor spec, parsed
  random effects).

- formula_details:

  Fixed-effect design (`X`, `rhs`, `contrasts`).

- subject_pars:

  Data frame of subject-level parameters. For a 1-RE fit (`k ~ 1`) the
  columns are `id, u_i, k`; for a phi-target 2-RE fit (`k + phi ~ 1`)
  they are `id, re_k, re_phi, k, phi, phi_latent` (`phi` floored at 0.1,
  `phi_latent` unfloored); for an s-target 2-RE fit (`k + s ~ 1`,
  GM/Rachlin) they are `id, re_k, re_s, k, s, s_latent` where `s` is
  soft-clamped toward `(0.05, 20)` and `s_latent = exp(log_s + re_s)`.

- guard_info:

  Guard / clamp / floor activity at the fitted values (see "Scales,
  guards and floors"): `n_rows`, `mu_guard_lower`, `mu_guard_upper`,
  plus `n_s_clamped_lower`/`n_s_clamped_upper` (`k + s ~ 1`) or
  `n_phi_floor` (`k + phi ~ 1`). `NULL` if the computation failed.

- loglik, AIC, BIC:

  Fit statistics.

- converged, se_available:

  Convergence / SE-availability flags. A non-converged fit raises a
  `beezdiscounting_convergence_warning` at fit time (regardless of
  `verbose`) and again from
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)
  and
  [`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md).

- multi_start_info:

  List recording the start selection: `n_starts`, `n_finite`,
  `n_converged`, `selected_start`, `tier` (1 = converged and passing the
  log-k sanity guard, 2 = passing the guard but not converged, 3 =
  neither), and `lower_nll_nonconverged` (`TRUE` when a non-converged
  start reached a lower NLL than the kept one).

- opt_warnings:

  Character vector of optimizer warnings.

- data:

  The single filtered model frame (id/x/y + retained design columns),
  row-aligned with the design matrix.

- data_all:

  The validated frame before complete-casing.

- coercion_info:

  Scale-coercion/clamping audit list.

## Scales, guards and floors

- **Random-effect scales.**
  [`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) reports
  the subject SD of `log k` on the natural-log scale (`Term = "k"` means
  log k); [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) convert it to the
  log10 scale and label it so.
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) returns
  the standardised deviate `u_i` for a `k ~ 1` fit but natural-log
  offsets (`re_k`, `re_phi` / `re_s`) for a two-random-effect fit.

- **Mean guard.** The fitted mean is held inside `[1e-6, 1 - 1e-6]`. For
  the exponential equation this binds once `k * delay` exceeds about
  13.8 (for `k = 0.01`, delays beyond roughly 1,400 days); such
  observations carry no information about `k`, so very steep discounters
  measured at long delays lose curvature. Mazur needs `k * delay` near
  `1e6` to reach the guard; Green-Myerson and Rachlin can reach it with
  a large `s` (e.g. Green-Myerson with `k = 0.1`, `s = 3` at 1,460
  days). `fit$guard_info$mu_guard_lower` / `$mu_guard_upper` count the
  positive-delay observations whose subject-level fitted mean falls
  outside the guard, and
  [`summary()`](https://rdrr.io/r/base/summary.html) adds a note when
  either is non-zero.

- **Shape `s`.** The reported `s` is the unclamped population value
  `exp(log_s)`. With `k + s ~ 1` each subject's effective `s` is
  soft-clamped into `(0.05, 20)`, and the `VarCorr()` SD of `s` is on
  the latent (pre-clamp) log scale. `subject_pars` holds both the
  effective `s` and the latent `s_latent = exp(log_s + re_s)`;
  `fit$guard_info$n_s_clamped_lower` / `$n_s_clamped_upper` count the
  subjects whose two values differ by more than 1% on the log scale
  (noted by [`summary()`](https://rdrr.io/r/base/summary.html)).

- **Precision floor.** SLT-beta precision is bounded below at
  `phi = 0.1`. For `k ~ 1` this is an optimizer bound that
  `tmb_control$lower` can relax; with `k + phi ~ 1` each subject's `phi`
  is floored at 0.1 inside the likelihood and cannot be relaxed.
  `subject_pars$phi_latent` is the unfloored value and
  `fit$guard_info$n_phi_floor` counts the subjects below the floor
  (noted by [`summary()`](https://rdrr.io/r/base/summary.html)).

These counts are computed at the fitted values and are reporting only:
they do not change any estimate.

## Two-parameter equations (Green-Myerson, Rachlin)

`k` and `s` trade off along a ridge, and with the few delays of a
typical titration task (about 7) the pair is much less stable than a
one-parameter `k`. The pair is particularly sensitive to how responses
near zero are recorded. In the SLT-beta log-density the response enters
through a `(mu * phi - 1) * log(y)` term, so a recorded 0 (treated as
about `1e-8`) carries far more leverage than a recorded 0.001; the
direction of the pull depends on `mu * phi`. In simulation with
identical true values (Green- Myerson, `log k = -4.61`, `s = 1.4`, 300
subjects, 6 delays), the same draws gave `(log k, s)` of `(-4.79, 1.58)`
as simulated, `(-5.06, 1.91)` after rounding `y` to 3 decimals, and
`(-4.37, 1.20)` after flooring `y` at 0.001; Mazur fits to the same data
moved by less than 0.03. These are sensitivity results for that design,
not a general bias. Report how the indifference points were recorded
(resolution, how zeros were coded, the count of `y` below 0.001), keep
the convention identical across groups being compared, and consider a
sensitivity refit under an alternative convention.

## References

Young, M. E. (2017). Discounting: A practical guide to multilevel
analysis of indifference data. *Journal of the Experimental Analysis of
Behavior, 108*(1), 97-112.
[doi:10.1002/jeab.265](https://doi.org/10.1002/jeab.265)

Kim, M., Koffarnus, M. N., & Franck, C. T. (2024). Thinking inside the
bounds: Improved error distributions for indifference point data
analysis and simulation via beta regression using common discounting
functions. *arXiv* preprint arXiv:2404.18000.

Kim, M., Kaplan, B. A., Koffarnus, M. N., & Franck, C. T. (2025).
Scale-location-truncated beta regression: Expanding beta regression to
accommodate 0 and 1. *arXiv* preprint arXiv:2509.13167.

## Examples

``` r
# \donttest{
# Small two-subject long-format indifference-point data frame.
dd <- data.frame(
  id = rep(c("s1", "s2"), each = 5),
  x  = rep(c(7, 30, 180, 365, 730), times = 2),
  y  = c(0.95, 0.80, 0.45, 0.30, 0.15,
         0.90, 0.70, 0.40, 0.25, 0.10)
)
fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb", verbose = 0)
exp(fit$model$coefficients[["beta_k"]])  # population k
#> [1] 0.008690582
# }
```
