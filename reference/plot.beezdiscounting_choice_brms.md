# Plot a Bayesian (brms) choice-based discounting model

Visualize a fitted
[`fit_dd_choice_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice_brms.md)
model. `type = "population"` (default) draws the posterior implied
discount curve with a credible band; `type = "individual"` adds
per-subject (posterior-median) curves; `type = "parameters"` shows the
posterior subject-`k` caterpillar; and `type = "calibration"` plots
observed choice proportion against fitted P(choose LL), binned by fitted
probability.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice_brms'
plot(
  x,
  type = c("population", "individual", "calibration", "parameters"),
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

  A `beezdiscounting_choice_brms` object.

- type:

  One of `"population"`, `"individual"`, `"calibration"`,
  `"parameters"`.

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

  Unused for the implied discount curve.

- ...:

  Unused.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
