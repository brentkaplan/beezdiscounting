# Glance at a beezdiscounting_brms model

The dd glance columns with `backend = "brms"` and
`logLik`/`AIC`/`BIC = NA_real_` by design, plus the Bayesian comparison
and diagnostic currency (`elpd_loo`, `p_loo`, `looic`, `rhat_max`,
`ess_bulk_min`, `num_divergent`).

## Usage

``` r
# S3 method for class 'beezdiscounting_brms'
glance(x, ...)
```

## Arguments

- x:

  A `beezdiscounting_brms` object.

- ...:

  Unused.

## Value

A one-row tibble.
