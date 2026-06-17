# Fit a trial-level SS-vs-LL choice model (binomial GLMM) via TMB

Fits a trial-level binary choice model from two complementary
perspectives. The structural model estimates the discount rate `k`
directly from choices via the scale-invariant relative-value comparison
`logit P(LL) = beta0 + gamma * ((ll/ss) * D(k, delay) - 1)`,
`k = exp(X beta_k + sigma_u u)`. The descriptive model (Young 2018)
forgoes a discount function and instead characterises choices via
separate magnitude and delay sensitivities with optional correlated
per-subject random slopes.

## Usage

``` r
fit_dd_choice(
  data,
  mode = c("structural", "descriptive"),
  id_var = "id",
  ss_var = "ss_amount",
  ll_var = "ll_amount",
  delay_var = "delay",
  choice_var = "choice",
  equation = c("mazur", "exponential"),
  intercept = FALSE,
  predictors = NULL,
  random_slopes = TRUE,
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  start_values = NULL,
  tmb_control = list(iter_max = 1000, eval_max = 2000),
  multi_start = TRUE,
  verbose = 1,
  ...
)
```

## Arguments

- data:

  Trial-level data frame (see the `*_var` args).

- mode:

  `"structural"` (default) estimates the discount rate `k` directly from
  choices via a discount function (Mazur or exponential); shares the IP
  family's `k`/emmeans contract. `"descriptive"` fits the Young (2018)
  correlated random-slope logistic model with separate magnitude and
  delay sensitivity predictors; the inferential targets are `VarCorr` /
  `ranef`, not emmeans.

- id_var, ss_var, ll_var, delay_var, choice_var:

  Column names.

- equation:

  `"mazur"` or `"exponential"` (structural mode only).

- intercept:

  Logical; include the choice-bias `beta0` (default `FALSE`). Structural
  only (the descriptive design carries its own intercept policy via
  `predictors`).

- predictors:

  Descriptive (`mode = "descriptive"`) only: a one-sided formula for the
  fixed-effect design on the logit of choosing the larger-later reward,
  or `NULL` for Young's (2018) two scale-invariant predictors
  (`~ 0 + log(ll_amount / ss_amount) + log(delay + 1)`). Ignored when
  `mode = "structural"`.

- random_slopes:

  Descriptive only: logical; `TRUE` (default) fits two correlated
  per-subject random slopes on Young's two predictors (`log(ll/ss)` and
  `log(delay + 1)`), giving a full bivariate random-effect covariance;
  `FALSE` fits a pooled fixed-effect logistic model with no random
  effects. Ignored when `mode = "structural"`.

- factors, factor_interaction, continuous_covariates:

  Between-subject design on `log k` (same semantics as
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)).
  Structural only.

- start_values, tmb_control, multi_start, verbose, ...:

  As in
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).

## Value

An object of class `beezdiscounting_choice`.

## Details

The **structural** model (`mode = "structural"`) parameterises choices
through a classical discount function (Mazur hyperbolic or exponential)
and shares the IP family's `k`/emmeans contract: `k` estimates and
estimated marginal means are accessible via
[`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)
and
[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md).
The **descriptive** model (`mode = "descriptive"`) follows Young (2018):
it regresses binary choice on `log(ll_amount / ss_amount)` (magnitude
sensitivity) and `log(delay + 1)` (delay sensitivity) with no assumed
discount function. When `random_slopes = TRUE` (default) each subject
receives correlated random slopes on these two predictors; the primary
inferential targets are the random-effect covariance (via
[`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html)) and
per-subject slopes (via
[`nlme::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html)),
not emmeans.

## See also

[`simulate_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/simulate_dd_choice.md)
for data generation;
[`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) and
[`nlme::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) for
descriptive-mode random-effect output;
[`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)
and
[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)
for structural-mode emmeans.

## Examples

``` r
# \donttest{
# Structural model: estimate a discount rate k from binary choices
sim_s <- simulate_dd_choice(n_subjects = 30, mode = "structural", seed = 1)
fit_s <- fit_dd_choice(sim_s, mode = "structural", equation = "mazur")
summary(fit_s)
#> 
#> Structural Choice Discounting Model Summary
#> ================================================== 
#> 
#> Call:
#> fit_dd_choice(data = sim_s, mode = "structural", equation = "mazur") 
#> 
#> Mode: structural 
#> Equation: mazur 
#> Backend: TMB_choice 
#> Convergence: Yes 
#> Subjects: 30  Observations: 300 
#> 
#> --- Fixed Effects (k) ---
#>           term estimate std.error statistic p.value
#>  k:(Intercept)   0.0217    0.0029   -29.156  <2e-16
#>          gamma   4.9956    0.7467    10.762  <2e-16
#> 
#> --- Variance Components ---
#>                Component Estimate Scale
#>  sigma_u (log10-k RE SD)   0.1339 log10
#> 
#> --- Fit Statistics ---
#> Log-likelihood: -113.6 
#> AIC: 233.19   BIC: 244.3 
#> 
#> Notes:
#>   * gamma is the choice-sensitivity (slope) parameter on the logit scale. 

# Descriptive (Young 2018) model: correlated per-subject magnitude/delay slopes
sim_d <- simulate_dd_choice(n_subjects = 30, mode = "descriptive", seed = 1)
fit_d <- fit_dd_choice(sim_d, mode = "descriptive")
summary(fit_d)
#> 
#> Descriptive (Young 2018) Choice Discounting Model Summary
#> ================================================== 
#> 
#> Call:
#> fit_dd_choice(data = sim_d, mode = "descriptive") 
#> 
#> Backend: TMB_choice   Convergence: Yes 
#> Subjects: 30  Observations: 300 
#> 
#> --- Fixed sensitivities (logit scale) ---
#>                      term estimate std.error statistic  p.value
#>  log(ll_amount/ss_amount)   1.3402    0.3576    3.7477 0.000178
#>            log(delay + 1)  -0.4633    0.3576   -1.2956 0.195110
#> 
#> --- Random-slope (co)variances ---
#>    Group                     Term   Variance    StdDev       Corr
#>  subject log(ll_amount/ss_amount) 0.11206266 0.3347576         NA
#>  subject           log(delay + 1) 0.04255008 0.2062767 -0.2443033
#> 
#> --- Fit Statistics ---
#> Log-likelihood: -172.67   AIC: 355.33   BIC: 373.85 
#> 
#> Notes:
#>   * theta are logit-scale (Young 2018) sensitivities; not exponentiated. 
nlme::VarCorr(fit_d)
#>     Group                     Term   Variance    StdDev       Corr
#> 1 subject log(ll_amount/ss_amount) 0.11206266 0.3347576         NA
#> 2 subject           log(delay + 1) 0.04255008 0.2062767 -0.2443033
# }
```
