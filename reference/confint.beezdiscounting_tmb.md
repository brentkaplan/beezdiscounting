# Confidence intervals for a TMB discounting model

Wald (Hessian-based) confidence intervals: `estimate +/- z * se` on the
internal (log-k) scale, then back-transformed to `report_space`.

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
confint(
  object,
  parm = NULL,
  level = 0.95,
  report_space = c("internal", "natural"),
  ...
)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- parm:

  Optional character vector for filtering. Accepts display names
  (`"k:(Intercept)"`, and `"s"` for the GM/Rachlin shape parameter)
  **or** raw optimizer names (`"beta_k"`, `"log_sigma_u"`, `"log_phi"`,
  `"log_sigma_e"`, `"log_s"`, and for a 2-RE fit `"log_sd_re"` /
  `"cor_re"`). `NULL` returns all coefficients.

- level:

  Confidence level (default `0.95`).

- report_space:

  `"internal"` (default; all coefficients on their estimation/log scale)
  or `"natural"` (exponentiate `beta_k` rows so the intercept is `k` at
  the reference level and non-intercept terms are multiplicative
  fold-changes; variance/aux params stay on their internal scale).

- ...:

  Unused.

## Value

A tibble with columns `term`, `estimate`, `conf.low`, `conf.high`,
`level`.

## Note

For a 2-RE fit the `cor_re` and `log_sd_re` rows are reported on their
internal (atanh / log) scales and are NOT back-transformed by
`report_space = "natural"` (only `beta_k`/`log_s` rows are); use
`VarCorr()` to obtain the correlation and natural-log SDs.

## Examples

``` r
# \donttest{
dd <- simulate_dd_ip(n_subjects = 20, seed = 1)
fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb",
                  random_effects = k ~ 1, verbose = 0)
confint(fit)
#> # A tibble: 3 × 5
#>   term          estimate conf.low conf.high level
#>   <chr>            <dbl>    <dbl>     <dbl> <dbl>
#> 1 k:(Intercept)   -4.53     -4.80    -4.27   0.95
#> 2 log_sigma_u     -0.632    -1.03    -0.235  0.95
#> 3 log_phi          2.52      2.29     2.76   0.95
confint(fit, parm = "k:(Intercept)", report_space = "natural")
#> # A tibble: 1 × 5
#>   term          estimate conf.low conf.high level
#>   <chr>            <dbl>    <dbl>     <dbl> <dbl>
#> 1 k:(Intercept)   0.0108  0.00826    0.0140  0.95
# }
```
