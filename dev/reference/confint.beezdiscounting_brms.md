# Credible intervals for a beezdiscounting_brms model

Equal-tailed quantile credible intervals on the report-space-transformed
draws. Columns match
[`confint.beezdiscounting_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/confint.beezdiscounting_tmb.md):
`term`, `estimate`, `conf.low`, `conf.high`, `level`.

## Usage

``` r
# S3 method for class 'beezdiscounting_brms'
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

  A `beezdiscounting_brms` object.

- parm:

  Optional terms (display, e.g. `"k:(Intercept)"`/`"s"`, or TMB
  coefficient names `"beta_k"`/`"log_s"`).

- level:

  Credible level (default 0.95).

- report_space:

  Reporting scale.

- ...:

  Unused.

## Value

A tibble.
