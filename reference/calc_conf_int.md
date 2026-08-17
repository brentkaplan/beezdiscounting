# Calculate Confidence Intervals for a Parameter

This function computes the lower and upper bounds of the confidence
interval for a parameter estimate, given its standard error, a specified
significance level, and the degrees of freedom from the model.

## Usage

``` r
calc_conf_int(estimate, std_error, model, alpha = 0.05)
```

## Arguments

- estimate:

  A numeric value representing the parameter estimate.

- std_error:

  A numeric value representing the standard error of the parameter
  estimate.

- model:

  A fitted model object that provides the residual degrees of freedom
  via [`df.residual()`](https://rdrr.io/r/stats/df.residual.html).

- alpha:

  A numeric value representing the significance level. Default is 0.05
  (95% confidence interval).

## Value

A numeric vector of length two:

- First element: Lower bound of the confidence interval.

- Second element: Upper bound of the confidence interval.

## Examples

``` r
# Example using a linear model
data <- data.frame(x = 1:10, y = c(2.3, 2.1, 3.7, 4.5, 5.1, 6.8, 7.3, 7.9, 9.2, 10.1))
lm_model <- lm(y ~ x, data = data)
calc_conf_int(estimate = 0.5, std_error = 0.1, model = lm_model, alpha = 0.05)
#> [1] 0.2693996 0.7306004
```
