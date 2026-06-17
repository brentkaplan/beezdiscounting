# Fitted values for a beezdiscounting_tmb fit

Fitted values for a beezdiscounting_tmb fit

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
fitted(object, level = c("subject", "population"), ...)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- level:

  `"subject"` (default) or `"population"`.

- ...:

  Unused.

## Value

Numeric vector of fitted indifference proportions, length
`nobs(object)`.

## Examples

``` r
# \donttest{
dd <- simulate_dd_ip(n_subjects = 20, seed = 1)
fit <- fit_dd_tmb(dd, verbose = 0)
head(fitted(fit))
#> [1] 0.95335370 0.82665521 0.44283742 0.28158867 0.16386570 0.08924494
# }
```
