# Plot a Bayesian (brms) indifference-point discounting model

Visualize a fitted
[`fit_dd_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_brms.md)
model. `type = "population"` (default) draws the posterior-median
discount curve with a credible band for the population conditional mean
(random effects set to zero – a credible band, not a
posterior-predictive band) over the observed (boundary-squeezed)
indifference points; `type = "individual"` adds per-subject median
curves; `type = "parameters"` shows the posterior subject-`k`
caterpillar; `type = "resid"` plots standardized residuals against
fitted values.

## Usage

``` r
# S3 method for class 'beezdiscounting_brms'
plot(
  x,
  type = c("population", "individual", "parameters", "resid"),
  ids = NULL,
  at = NULL,
  n_points = 200,
  x_trans = c("log10", "linear"),
  show_observed = TRUE,
  ...
)
```

## Arguments

- x:

  A `beezdiscounting_brms` object.

- type:

  One of `"population"`, `"individual"`, `"parameters"`, `"resid"`.

- ids:

  Optional subset of subject ids for `type = "individual"`.

- at:

  Optional named list conditioning the population curve on factor levels
  / covariate values.

- n_points:

  Number of delay points in the curve grid.

- x_trans:

  Delay-axis scale: `"log10"` (default) or `"linear"`.

- show_observed:

  Overlay the observed indifference points.

- ...:

  Unused.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
