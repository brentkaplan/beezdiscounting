# Residuals for a beezdiscounting_tmb fit

Residuals for a beezdiscounting_tmb fit

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
residuals(
  object,
  type = c("response", "pearson"),
  level = c("subject", "population"),
  ...
)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- type:

  `"response"` (default; raw `y - fitted`) or `"pearson"` (divided by
  the per-row response SD; equivalent to `.std_resid` in
  [`augment()`](https://generics.r-lib.org/reference/augment.html)).

- level:

  `"subject"` (default) or `"population"`.

- ...:

  Unused.

## Value

Numeric vector of residuals, length `nobs(object)`.

## Examples

``` r
# \donttest{
dd <- simulate_dd_ip(n_subjects = 20, seed = 1)
fit <- fit_dd_tmb(dd, verbose = 0)
head(residuals(fit))
#> [1]  0.044642112  0.064130983  0.127092182  0.008255303 -0.011513496
#> [6]  0.056085650
head(residuals(fit, type = "pearson"))
#> [1]  0.77723175  0.62200194  0.93939263  0.06738764 -0.11420013  0.72227170
# }
```
