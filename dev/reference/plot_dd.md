# Plot Delay-Discounting Model

This function generates a plot of the delay-discounting data and the
fitted model.

## Usage

``` r
plot_dd(
  fit_dd_object,
  xlabel = "Delay",
  ylabel = "Indifference Point",
  title = "",
  logx = TRUE
)
```

## Arguments

- fit_dd_object:

  A fitted delay-discounting model object of class `"fit_dd"`, created
  by the
  [`fit_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd.md)
  function.

- xlabel:

  A character string specifying the label for the x-axis. Default is
  `"Delay"`.

- ylabel:

  A character string specifying the label for the y-axis. Default is
  `"Indifference Point"`.

- title:

  A character string specifying the plot title. Default is `""`.

- logx:

  Logical. If `TRUE`, the x-axis is log-transformed. Default is `TRUE`.

## Value

A ggplot object representing the fitted model and data.

## Examples

``` r
data <- data.frame(
  id = rep(1:2, each = 6),
  x = rep(c(1, 7, 30, 90, 180, 365), 2),
  y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
)
fit <- fit_dd(data, equation = "mazur", method = "mean")
plot_dd(fit)
```
