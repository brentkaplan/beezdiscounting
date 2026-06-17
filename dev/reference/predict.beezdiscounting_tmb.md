# Predict from a TMB mixed-effects discounting model

Predict from a TMB mixed-effects discounting model

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
predict(
  object,
  newdata = NULL,
  type = c("response", "parameters"),
  level = "subject",
  ...
)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- newdata:

  Optional data frame. `NULL` uses the fitting data. A supplied
  `newdata` must use the package's canonical column names (`id`, `x`,
  `y`), regardless of the `id_var`/`x_var`/`y_var` you gave at fit time.

- type:

  `"response"` (fitted indifference proportions on `(0,1)`) or
  `"parameters"` (the per-subject parameter tibble).

- level:

  For `type = "response"`: `"subject"` (default; conditions on each
  subject's estimated random intercept, requires the id column) and/or
  `"population"` (random effects set to zero - the population-mean
  curve; no id column needed). Pass `c("population", "subject")` for
  both columns side-by-side. A numeric nlme-style level is rejected with
  an error. For an s-target 2-RE fit (`k + s ~ 1`), the subject level
  uses each subject's estimated `s_i`; the population level (random
  effects zero) uses the soft-clamped population
  `s = .dd_soft_clamp_s_log(log_s)`, matching the kernel near the
  `(0.05, 20)` bounds (equal to `exp(log_s)` in the interior).

- ...:

  Unused.

## Value

- `type = "parameters"`: the subject-parameter tibble (one row per
  subject).

- `type = "response"`, `level = "subject"` only: `newdata` as a tibble
  plus a `.fitted` column.

- `type = "response"`, `"population"` requested: `newdata` plus
  `predict.fixed` (and `predict.id` when `"subject"` is also requested),
  matching the `nlme::predict.lme(level = 0:1)` column-name convention
  so nlme-based plotting code runs unchanged.

## Examples

``` r
# \donttest{
dd <- simulate_dd_ip(n_subjects = 20, seed = 1)
fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb",
                  random_effects = k ~ 1, verbose = 0)

# Subject-conditional fitted values
head(predict(fit, type = "response"))
#> # A tibble: 6 × 4
#>   id        x     y .fitted
#>   <chr> <dbl> <dbl>   <dbl>
#> 1 1         7 0.998  0.953 
#> 2 1        30 0.891  0.827 
#> 3 1       180 0.570  0.443 
#> 4 1       365 0.290  0.282 
#> 5 1       730 0.152  0.164 
#> 6 1      1460 0.145  0.0892

# Population-mean curve at specific delays (no id column needed)
nd <- data.frame(x = c(7, 30, 180, 365))
predict(fit, newdata = nd, level = "population")
#> # A tibble: 4 × 2
#>       x predict.fixed
#>   <dbl>         <dbl>
#> 1     7         0.930
#> 2    30         0.756
#> 3   180         0.341
#> 4   365         0.203

# Per-subject parameters
head(predict(fit, type = "parameters"))
#> # A tibble: 6 × 3
#>   id        u_i       k
#>   <chr>   <dbl>   <dbl>
#> 1 1     -0.811  0.00699
#> 2 10    -0.709  0.00738
#> 3 11     1.33   0.0218 
#> 4 12     0.0467 0.0110 
#> 5 13    -0.543  0.00806
#> 6 14    -2.30   0.00317
# }
```
