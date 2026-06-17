# Extract subject-level random effects from a choice model

Structural: `id`, the standardized random-intercept deviate `u_i` (such
that `log k_i = X beta + sigma_u * u_i`), and the resolved per-subject
discount rate `k`. Descriptive (Young 2018): `id` plus the per-subject
magnitude and delay random slopes (`b_mag`, `b_delay`); pooled
(`random_slopes = FALSE`) returns an id-only frame.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
ranef(object, ...)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

- ...:

  Unused.

## Value

Data frame: structural `id`/`u_i`/`k`; descriptive `id`/`b_mag`/
`b_delay` (or id-only when pooled).
