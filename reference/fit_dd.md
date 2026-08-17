# Fit Delay-Discounting Model

This function fits a delay-discounting model to the given dataset using
the specified equation and method.

## Usage

``` r
fit_dd(dat, equation, method)
```

## Arguments

- dat:

  A data frame containing delay (`x`) and indifference point (`y`) data.
  For "two stage" methods, the data must include an `id` column to
  identify participants.

- equation:

  A character string specifying the delay-discounting equation to use.
  Options include:

  - `"mazur"` or `"hyperbolic"`: Hyperbolic delay-discounting model (\\y
    = 1 / (1 + k \cdot x)\\).

  - `"exponential"`: Exponential delay-discounting model (\\y = \exp(-k
    \cdot x)\\).

- method:

  A character string specifying the method for fitting the model.
  Options include:

  - `"pooled"` or `"agg"`: Fits the model using pooled data.

  - `"mean"`: Fits the model using the mean of indifference points at
    each delay.

  - `"ts"` or `"two stage"`: Fits the model separately for each
    participant (requires an `id` column in `dat`).

## Value

A list object of class `"fit_dd"`, containing:

- The fitted model(s).

- The original dataset (`dat`).

- The specified method (`method`).

## Examples

``` r
data <- data.frame(
  id = rep(1:2, each = 6),
  x = rep(c(1, 7, 30, 90, 180, 365), 2),
  y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
)
fit_dd(data, equation = "mazur", method = "two stage")
#> [[1]]
#> [[1]]$`1`
#> [[1]]$`1`$result
#> Nonlinear regression model
#>   model: y ~ 1/(1 + k * x)
#>    data: data
#>      k 
#> 0.1008 
#>  residual sum-of-squares: 0.02315
#> 
#> Number of iterations to convergence: 10 
#> Achieved convergence tolerance: 1.49e-08
#> 
#> [[1]]$`1`$error
#> NULL
#> 
#> 
#> [[1]]$`2`
#> [[1]]$`2`$result
#> Nonlinear regression model
#>   model: y ~ 1/(1 + k * x)
#>    data: data
#>       k 
#> 0.08777 
#>  residual sum-of-squares: 0.01866
#> 
#> Number of iterations to convergence: 9 
#> Achieved convergence tolerance: 1.49e-08
#> 
#> [[1]]$`2`$error
#> NULL
#> 
#> 
#> 
#> [[2]]
#>    id   x    y
#> 1   1   1 0.90
#> 2   1   7 0.50
#> 3   1  30 0.30
#> 4   1  90 0.20
#> 5   1 180 0.10
#> 6   1 365 0.05
#> 7   2   1 0.85
#> 8   2   7 0.55
#> 9   2  30 0.35
#> 10  2  90 0.15
#> 11  2 180 0.10
#> 12  2 365 0.05
#> 
#> [[3]]
#> [1] "two stage"
#> 
#> attr(,"class")
#> [1] "fit_dd" "list"  
```
