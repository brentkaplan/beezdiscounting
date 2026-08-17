# Extract Results from Delay-Discounting Model

This function extracts model parameter estimates, fit statistics, and
confidence intervals from a fitted delay-discounting model.

## Usage

``` r
results_dd(fit_dd_object)
```

## Arguments

- fit_dd_object:

  A fitted delay-discounting model object of class `"fit_dd"`, created
  by the
  [`fit_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd.md)
  function.

## Value

A tibble containing the following columns:

- `id`: The participant or group ID (if applicable).

- `term`: The model parameter (e.g., `k`).

- `estimate`: The estimated value of the parameter.

- `std.error`: The standard error of the parameter estimate.

- `statistic`: The t-statistic for the parameter estimate.

- `p.value`: The p-value for the parameter estimate.

- `conf_low`: The lower bound of the 95% confidence interval.

- `conf_high`: The upper bound of the 95% confidence interval.

- `R2`: The coefficient of determination (\\R^2\\).

## Examples

``` r
data <- data.frame(
  id = rep(1:2, each = 6),
  x = rep(c(1, 7, 30, 90, 180, 365), 2),
  y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
)
fit <- fit_dd(data, equation = "mazur", method = "two stage")
results_dd(fit)
#> # A tibble: 2 × 22
#>   method  id    term  estimate std.error statistic p.value  sigma isConv  finTol
#>   <chr>   <chr> <chr>    <dbl>     <dbl>     <dbl>   <dbl>  <dbl> <lgl>    <dbl>
#> 1 two st… 1     k       0.101     0.0205      4.91 0.00443 0.0680 TRUE   1.49e-8
#> 2 two st… 2     k       0.0878    0.0158      5.55 0.00262 0.0611 TRUE   1.49e-8
#> # ℹ 12 more variables: logLik <dbl>, AIC <dbl>, BIC <dbl>, deviance <dbl>,
#> #   df.residual <int>, nobs <int>, R2 <dbl>, auc_regular <dbl>,
#> #   auc_log10 <dbl>, auc_ord <dbl>, conf_low <dbl>, conf_high <dbl>
```
