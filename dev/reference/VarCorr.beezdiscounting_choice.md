# Random-effect (co)variances for a beezdiscounting_choice fit

Descriptive (Young 2018): the 2-RE slope SDs and their correlation (the
per-subject magnitude/delay sensitivities); pooled
(`random_slopes = FALSE`) returns a 0-row frame. Structural: the log-k
random-intercept SD.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
VarCorr(x, sigma = 1, ...)
```

## Arguments

- x:

  A `beezdiscounting_choice` object.

- sigma:

  Ignored (present for
  [`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) generic
  compatibility).

- ...:

  Unused.

## Value

A data frame of variance components. Descriptive: `Group`, `Term`,
`Variance`, `StdDev`, `Corr` (one row per random slope; 0 rows when
pooled). Structural: `Group`, `Term`, `StdDev`, `Scale`.
