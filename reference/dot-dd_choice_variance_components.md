# Variance components for a beezdiscounting_choice fit

Reports the random-intercept SD on the log10-k scale (divides the
natural-log SD by `log(10)` for comparability with
[`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html)).

## Usage

``` r
.dd_choice_variance_components(object)
```

## Arguments

- object:

  A `beezdiscounting_choice` fit.

## Value

A data frame with columns `Component`, `Estimate`, `Scale`.
