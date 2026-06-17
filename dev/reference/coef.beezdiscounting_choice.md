# Extract coefficients from a structural choice model

Returns the optimizer's flat named parameter vector: `beta_k` (one entry
per fixed-effect design column, on the log-k scale), `log_sigma_u`,
`log_gamma`, and `beta0` (only when `intercept = TRUE`). This is the
numeric escape hatch consumed by tooling such as `car::deltaMethod`.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
coef(object, ...)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

- ...:

  Unused.

## Value

Named numeric vector.
