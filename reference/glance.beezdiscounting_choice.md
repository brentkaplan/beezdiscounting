# Glance at a beezdiscounting_choice model

Returns a one-row tibble of model-level summary statistics suitable for
use with
[`broom::glance()`](https://generics.r-lib.org/reference/glance.html)
workflows and multi-model comparison tables.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
glance(x, ...)
```

## Arguments

- x:

  A `beezdiscounting_choice` object.

- ...:

  Unused.

## Value

A one-row tibble with columns `model_class`, `backend`, `mode`,
`equation`, `nobs`, `n_subjects`, `n_random_effects`, `converged`,
`logLik`, `AIC`, `BIC`.
