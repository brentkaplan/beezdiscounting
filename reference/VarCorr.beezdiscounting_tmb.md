# Random-effect covariance for a TMB discounting model

For a 2-RE fit returns the 2x2 covariance with `StdDev` and the joint
correlation (a structural `0` for `pdDiag`; the first row's `Corr` is
`NA` by convention). The `Term` column reflects the fitted pair:
`(k, phi)` for a phi-target fit (`k + phi ~ 1`) or `(k, s)` for an
s-target fit (`k + s ~ 1`, GM/Rachlin only). For a 1-RE fit (`k ~ 1`)
returns the single random-intercept SD on the log-k scale.

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
VarCorr(x, sigma = 1, ...)
```

## Arguments

- x:

  A `beezdiscounting_tmb` fit.

- sigma:

  Ignored (present for the
  [`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html)
  generic).

- ...:

  Unused.

## Value

A data frame with `Group`, `Term`, `Variance`, `StdDev`, and (2-RE)
`Corr`.
