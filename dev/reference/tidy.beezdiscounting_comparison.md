# Tidy a discounting comparison into a flat contrasts frame

[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
method for `beezdiscounting_comparison` objects (returned by
[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)).
Produces a flat long tibble whose **column names and order** match the
beezdemand `beezdemand_comparison` tidier, enabling downstream consumers
to bind rows from both backends into a single frame. Note that the
`contrast` label dialect differs: beezdiscounting emits native
`factor=level` labels (e.g. `"condition=C1 - condition=C2"`), while
beezdemand emits bare level labels (e.g. `"C1 - C2"`). Fully uniform
emmeans-style contrast labels across backends are a future cross-backend
item. The nested object keeps the native dialect (see
[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)).

## Usage

``` r
# S3 method for class 'beezdiscounting_comparison'
tidy(x, exponentiate = FALSE, ...)
```

## Arguments

- x:

  A `beezdiscounting_comparison` object.

- exponentiate:

  Logical. If `TRUE`, return base-invariant ratios
  (`estimate = 10^estimate`, CIs back-transformed); `std.error` becomes
  `NA` following broom's convention for exponentiated fits. Default
  `FALSE` (log10-scale contrasts).

- ...:

  Unused.

## Value

A tibble with columns `param`, `contrast`, `estimate`, `std.error`,
`statistic`, `df`, `conf.low`, `conf.high`, `p.value` (with leading
by-column(s) inserted before `param` when `contrast_by` is active).
Estimates and CIs are on the log10 scale (or ratios when
`exponentiate = TRUE`). `statistic` is an asymptotic *z* (`df = Inf`).
On the brms backend, `statistic`/`df`/`p.value` are `NA` (posterior
summaries), the intervals are quantile credible intervals, and an
additional `post.prob` column reports the posterior probability of
direction.
