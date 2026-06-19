# Plot group differences in discount rate

Forest plot of the pairwise (or treatment-vs-control) contrasts from
[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md).
`type = "ratio"` (default) shows each k ratio on a log axis with a
reference line at 1; `type = "difference"` shows the difference in log10
k with a reference line at 0. A contrast is flagged when its interval
excludes the null – a backend-agnostic encoding that works for both the
frequentist (TMB) and Bayesian (brms) backends.

## Usage

``` r
# S3 method for class 'beezdiscounting_comparison'
plot(x, type = c("ratio", "difference"), ...)
```

## Arguments

- x:

  A `beezdiscounting_comparison` object.

- type:

  `"ratio"` (default) or `"difference"`.

- ...:

  Unused.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
# \donttest{
sim <- simulate_dd_ip(n_subjects = 30, n_conditions = 2, delta_k = c(0, 0.5), seed = 1)
fit <- fit_dd_tmb(sim, equation = "mazur", factors = "condition")
#> Fitting TMB mixed-effects discounting model (mazur, sltb)...
#>   Subjects: 30, Observations: 210
#>   Multi-start: best NLL = -339.65 (start set 1 of 3)
#>   Converged (NLL = -339.65). Done.
cmp <- get_dd_comparisons(fit)
plot(cmp)

# }
```
