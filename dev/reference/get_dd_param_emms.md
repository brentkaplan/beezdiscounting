# Estimated marginal means of the discount rate `k`

Computes estimated marginal means (EMMs) of the discount rate `k` from a
fitted `beezdiscounting_tmb` model (or a structural
`beezdiscounting_choice` fit, which shares the `log k` design). EMMs are
computed on the `log k` scale (linear in the fixed-effect coefficients)
using the averaging-matrix reference grid, then back-transformed with
[`exp()`](https://rdrr.io/r/base/Log.html) so that `k = exp(k_log)`.
Standard errors use the `beta_k` block of
[`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html)'s
fixed-effect covariance; intervals are Wald on the log scale and
exponentiated. Bayesian fits (`beezdiscounting_brms`,
`beezdiscounting_choice_brms`) use the same reference grid with
draws-based summaries: posterior medians and equal-tailed quantile
intervals, no delta method.

## Usage

``` r
get_dd_param_emms(fit, factors_in_emm = NULL, at = NULL, ci_level = 0.95, ...)
```

## Arguments

- fit:

  A `beezdiscounting_tmb`, structural `beezdiscounting_choice`,
  `beezdiscounting_brms`, or `beezdiscounting_choice_brms` object.

- factors_in_emm:

  Character vector of factors to retain in the EMM reference grid. A
  strict subset marginalizes the omitted factors with equal weights
  across the full crossing of their levels (emmeans' default
  `weights = "equal"`); `NULL` (default) retains all fitted factors;
  `character(0)` marginalizes everything to a single grand-mean cell.

- at:

  Named list specifying factor levels and/or continuous-covariate values
  for conditional EMMs (one numeric per covariate; multiple values warn
  and use the first). `at` on an omitted factor restricts the level set
  averaged over.

- ci_level:

  Numeric confidence level for intervals (default 0.95).

- ...:

  Additional arguments (currently unused).

## Value

A tibble with columns `level`, `k`, `k_log`, `std.error`, `conf.low`,
`conf.high`. `k_log` is the marginal mean on the `log k` scale;
`k = exp(k_log)`; `std.error` is the SE of `k_log`; the intervals are on
the `k` (natural) scale.

## Examples

``` r
# \donttest{
sim <- simulate_dd_ip(
  n_subjects = 30, n_conditions = 2, delta_k = c(0, log(3)), seed = 1
)
fit <- fit_dd_tmb(sim, factors = "condition", verbose = 0)
get_dd_param_emms(fit)
#> # A tibble: 2 × 6
#>   level             k k_log std.error conf.low conf.high
#>   <chr>         <dbl> <dbl>     <dbl>    <dbl>     <dbl>
#> 1 condition=C1 0.0112 -4.49     0.152  0.00834    0.0151
#> 2 condition=C2 0.0276 -3.59     0.154  0.0204     0.0373
# }
```
