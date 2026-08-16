# Plot questionnaire scores (internal implementation)

Shared boxplot implementation for
[`plot.score_mcq27_output()`](https://brentkaplan.github.io/beezdiscounting/reference/plot.score_mcq27_output.md),
[`plot.score_mcq_output()`](https://brentkaplan.github.io/beezdiscounting/reference/plot.score_mcq_output.md),
and
[`plot.score_pdq_output()`](https://brentkaplan.github.io/beezdiscounting/reference/plot.score_pdq_output.md).

## Usage

``` r
.plot_score_mcq(
  x,
  xlab = "Metric",
  alpha = 0.3,
  target_levels = c("small_k", "medium_k", "large_k", "geomean_k", "overall_k"),
  param = "k"
)
```

## Arguments

- x:

  A data frame returned by
  [`score_mcq()`](https://brentkaplan.github.io/beezdiscounting/reference/score_mcq.md),
  [`score_mcq27()`](https://brentkaplan.github.io/beezdiscounting/reference/score_mcq27.md),
  or
  [`score_pdq()`](https://brentkaplan.github.io/beezdiscounting/reference/score_pdq.md).

- xlab:

  Label for the x-axis.

- alpha:

  Transparency of the points in the plot.

- target_levels:

  Metric column names, in display order.

- param:

  Parameter letter used in column suffixes and the y label ("k" or "h").

## Value

A ggplot object showing the boxplot of scores.
