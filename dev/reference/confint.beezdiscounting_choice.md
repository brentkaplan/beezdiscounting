# Confidence intervals for a structural choice discounting model

Wald (Hessian-based) confidence intervals: `estimate +/- z * se` on the
internal scale, then back-transformed to `report_space`.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
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

  A `beezdiscounting_choice` object.

- parm:

  Optional character vector for filtering. Accepts display names
  (`"k:(Intercept)"`, `"gamma"`) **or** raw optimizer names (`"beta_k"`,
  `"log_gamma"`, `"log_sigma_u"`, `"beta0"`). `NULL` returns all
  coefficients.

- level:

  Confidence level (default `0.95`).

- report_space:

  `"internal"` (default; all coefficients on their estimation scale) or
  `"natural"` (exponentiate the `beta_k` rows so the intercept is `k` at
  the reference level, and exponentiate `log_gamma` so the row is
  `gamma`; `log_sigma_u` and `beta0` stay on their internal scales).

- ...:

  Unused.

## Value

A tibble with columns `term`, `estimate`, `conf.low`, `conf.high`,
`level`.
