# Extract subject-level random effects from a TMB discounting model

Extract subject-level random effects from a TMB discounting model

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
ranef(object, ...)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- ...:

  Unused.

## Value

Data frame keyed by `id`. For a 1-RE fit (`k ~ 1`): the standardized
random-intercept deviate `u_i` (such that
`log k_i = X beta + sigma_u * u_i`) and the resolved per-subject
discount rate `k` (no `phi` or `s` column - both are population-level).
For a phi-target 2-RE fit (`k + phi ~ 1`): `(re_k, re_phi)` offsets plus
per-subject `k` and `phi`. For an s-target 2-RE fit (`k + s ~ 1`,
GM/Rachlin only): `(re_k, re_s)` offsets plus per-subject `k` and `s`
(soft-clamped toward `(0.05, 20)`).
