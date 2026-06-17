# Plot Proportion of SIR/SS Choices by k Value

This function creates a plot of the proportion of SIR/SS choices by k
value using the output of the `prop_ss` function.

## Usage

``` r
# S3 method for class 'prop_ss_output'
plot(
  x,
  ...,
  pt_shape = 21,
  pt_fill = "white",
  pt_size = 3,
  title = "Proportion of SIR/SS choices by k value",
  xlab = "k value rank",
  ylab = "Proportion of SS choices"
)
```

## Arguments

- x:

  Output from the `prop_ss` function

- ...:

  Additional arguments passed to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

- pt_shape:

  Shape of the points in the plot. Default is 21.

- pt_fill:

  Fill color of the points in the plot. Default is "white".

- pt_size:

  Size of the points in the plot. Default is 3.

- title:

  Title of the plot. Default is "Proportion of SIR/SS choices by k
  value".

- xlab:

  Label for the x-axis. Default is "k value rank".

- ylab:

  Label for the y-axis. Default is "Proportion of SS choices".

## Value

A ggplot object.

## Examples

``` r
plot(prop_ss(mcq27))
```
