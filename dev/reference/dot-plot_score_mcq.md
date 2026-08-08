# Plot MCQ scores (internal implementation)

Shared boxplot implementation for
[`plot.score_mcq27_output()`](https://brentkaplan.github.io/beezdiscounting/reference/plot.score_mcq27_output.md)
and
[`plot.score_mcq_output()`](https://brentkaplan.github.io/beezdiscounting/reference/plot.score_mcq_output.md).

## Usage

``` r
.plot_score_mcq(x, xlab = "Metric", alpha = 0.3)
```

## Arguments

- x:

  A data frame returned by
  [`score_mcq()`](https://brentkaplan.github.io/beezdiscounting/reference/score_mcq.md),
  [`score_mcq27()`](https://brentkaplan.github.io/beezdiscounting/reference/score_mcq27.md).

- xlab:

  Label for the x-axis.

- alpha:

  Transparency of the points in the plot.

## Value

A ggplot object showing the boxplot of MCQ scores.
