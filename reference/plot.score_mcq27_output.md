# Plot MCQ-27 Scores

Boxplots the MCQ-27 scores for different metrics (small_k, medium_k,
large_k, geomean_k, overall_k). Log-transformed k values
(`trans = "log"`/`"ln"`) are detected and the y-axis label adjusted
accordingly.

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
