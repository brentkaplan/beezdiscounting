# Augment a beezdiscounting_tmb model

Returns the data used for fitting (or `newdata`) as a tibble with three
diagnostic columns appended:

- `.fitted`: subject-conditional fitted indifference proportion (clamped
  to `(0, 1)`).

- `.resid`: raw residual `y - .fitted` on the response scale.

- `.std_resid`: Pearson (standardized) residual - `.resid` divided by
  the per-row response SD. For `family = "gaussian"` the SD is the
  constant `sigma_e`; for `family = "sltb"` it is the delta-method SLT
  SD `sqrt(mu * (1 - mu) / (phi + 1))`.

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
augment(x, newdata = NULL, ...)
```

## Arguments

- x:

  A `beezdiscounting_tmb` object.

- newdata:

  Optional data frame. `NULL` uses the fitting data. A supplied
  `newdata` must use the package's canonical column names (`id`, `x`,
  `y`), regardless of the `id_var`/`x_var`/`y_var` you gave at fit time.

- ...:

  Unused.

## Value

A tibble with the same rows as the data plus `.fitted`, `.resid`, and
`.std_resid`.

## Examples

``` r
# \donttest{
dd <- simulate_dd_ip(n_subjects = 20, seed = 1)
fit <- fit_dd_tmb(dd, verbose = 0)
head(augment(fit))
#> # A tibble: 6 × 6
#>   id        x     y .fitted   .resid .std_resid
#>   <chr> <dbl> <dbl>   <dbl>    <dbl>      <dbl>
#> 1 1         7 0.998  0.953   0.0446      0.777 
#> 2 1        30 0.891  0.827   0.0641      0.622 
#> 3 1       180 0.570  0.443   0.127       0.939 
#> 4 1       365 0.290  0.282   0.00826     0.0674
#> 5 1       730 0.152  0.164  -0.0115     -0.114 
#> 6 1      1460 0.145  0.0892  0.0561      0.722 
# }
```
