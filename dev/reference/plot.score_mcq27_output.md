# Plot MCQ-27 Scores

This function creates a plot of the MCQ-27 scores for different metrics
(small_k, medium_k, large_k, geomean_k, overall_k). The function handles
different logarithmic transformations of the k-values and adjusts the
y-axis label accordingly.

## Usage

``` r
# S3 method for class 'score_mcq27_output'
plot(x, ..., xlab = "Metric", alpha = 0.3)
```

## Arguments

- x:

  A data frame returned by the `score_mcq27` function.

- ...:

  Additional arguments passed to methods.

- xlab:

  Label for the x-axis. Default is "Metric".

- alpha:

  Transparency of the points in the plot. Default is 0.3.

## Value

A ggplot object showing the boxplot of MCQ-27 scores.

## Examples

``` r
plot(score_mcq27(mcq27))
```
