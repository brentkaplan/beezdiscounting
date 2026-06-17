# Extract coefficients from a TMB discounting model

Returns the optimizer's flat named parameter vector: `beta_k` (one entry
per fixed-effect design column, on the log-k scale), `log_sigma_u`, and
the auxiliary parameter (`log_phi` for `family = "sltb"`, `log_sigma_e`
for `family = "gaussian"`). This is the numeric escape hatch consumed by
tooling such as `car::deltaMethod`.

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
coef(object, ...)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- ...:

  Unused.

## Value

Named numeric vector.
