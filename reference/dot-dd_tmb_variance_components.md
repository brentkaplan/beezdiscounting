# Random-effect and auxiliary variance components for a beezdiscounting_tmb fit

Reports the random-intercept SD on the log10-k scale (divides the
natural-log SD by `log(10)` for comparability with
[`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html)) and the
auxiliary parameter back-transformed to its natural scale: `phi` (sltb)
or `sigma_e` (gaussian).

## Usage

``` r
.dd_tmb_variance_components(object)
```

## Arguments

- object:

  A `beezdiscounting_tmb` fit.

## Value

A data frame with columns `Component`, `Estimate`, `Scale`.
