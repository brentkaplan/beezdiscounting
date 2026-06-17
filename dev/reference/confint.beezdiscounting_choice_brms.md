# Credible intervals for a beezdiscounting_choice_brms model

Quantile credible intervals on the report-space-transformed draws;
columns `term`, `estimate`, `conf.low`, `conf.high`, `level`. `beta0`
stays on the identity scale across report spaces, matching
[`tidy.beezdiscounting_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_choice.md).

## Usage

``` r
# S3 method for class 'beezdiscounting_choice_brms'
confint(
  object,
  parm = NULL,
  level = 0.95,
  report_space = c("natural", "log10", "internal", "log"),
  ...
)
```

## Arguments

- object:

  A `beezdiscounting_choice_brms` object.

- parm:

  Optional terms: display names (`"k:(Intercept)"`, `"k:<level>"`,
  `"gamma"`, `"beta0"`) or TMB coefficient names (`"beta_k"`, which
  selects every log-k design coefficient, or `"log_gamma"`).

- level:

  Credible level (default 0.95).

- report_space:

  Reporting scale.

- ...:

  Unused.

## Value

A tibble.
