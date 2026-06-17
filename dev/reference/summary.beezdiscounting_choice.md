# Summarize a structural choice discounting fit

Summarize a structural choice discounting fit

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
summary(object, report_space = c("natural", "log10", "internal", "log"), ...)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

- report_space:

  Scale for the fixed-effect (`beta_k`) estimates and standard errors in
  the coefficient table: `"natural"` (default; `k` via
  [`exp()`](https://rdrr.io/r/base/Log.html)), `"log10"`, `"internal"`
  or `"log"`. `statistic` and `p.value` are always computed on the
  estimation scale regardless of `report_space`.

- ...:

  Unused.

## Value

An object of class `"summary.beezdiscounting_choice"` with components:
`call`, `model_class`, `backend`, `mode`, `equation`, `coefficients`,
`variance_components`, `n_subjects`, `nobs`, `converged`, `logLik`,
`AIC`, `BIC`, `notes`.
