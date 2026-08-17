# Tidy a beezdiscounting_choice_brms model

The 8-column dd contract; `gamma` (choice sensitivity) and `k` rows are
report-space transformed per draw; the optional `beta0` bias stays on
the logit scale.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice_brms'
tidy(
  x,
  effects = c("fixed", "ran_pars"),
  report_space = c("natural", "log10", "internal", "log"),
  ...
)
```

## Arguments

- x:

  A `beezdiscounting_choice_brms` object.

- effects:

  `"fixed"`, `"ran_pars"`, or both.

- report_space:

  Reporting scale.

- ...:

  Unused.

## Value

A tibble.
