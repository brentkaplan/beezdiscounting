# Calculate R-Squared for a Model

This function calculates the coefficient of determination (\\R^2\\) for
a given model by comparing the sum of squared errors (SSE) to the total
sum of squares (SST).

## Usage

``` r
calc_r2(model)
```

## Arguments

- model:

  A fitted model object. The model must have
  [`resid()`](https://rdrr.io/r/stats/residuals.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) methods to
  extract residuals and fitted values.

## Value

A numeric value representing the \\R^2\\ value of the model. Returns
`NA` if the model is `NULL`.

## Examples

``` r
# Example using a simple linear model
data <- data.frame(x = 1:10, y = c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11))
lm_model <- lm(y ~ x, data = data)
calc_r2(lm_model)
#> [1] 0.9927686
```
