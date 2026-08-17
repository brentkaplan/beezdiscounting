# Tidy a beezdiscounting_brms model

The exact 8-column dd coefficient contract (`term`, `estimate`,
`std.error`, `statistic`, `p.value`, `component`, `estimate_scale`,
`term_display`). Estimates are posterior medians and `std.error`
posterior SDs of the report-space-transformed draws (exact; no delta
method); `statistic`/`p.value` are `NA` – use
[`confint()`](https://rdrr.io/r/stats/confint.html) for interval
summaries. The shape row (`s`, two-parameter equations) carries
`component = "shape"`; variance rows (`effects = "ran_pars"`) mirror the
TMB reporting convention (log10-scale k RE SD; natural phi/sigma).

## Usage

``` r
# S3 method for class 'beezdiscounting_brms'
tidy(
  x,
  effects = c("fixed", "ran_pars"),
  report_space = c("natural", "log10", "internal", "log"),
  ...
)
```

## Arguments

- x:

  A `beezdiscounting_brms` object.

- effects:

  `"fixed"`, `"ran_pars"`, or both (default).

- report_space:

  `"natural"` (default), `"log10"`, `"internal"`, or `"log"`.

- ...:

  Unused.

## Value

A tibble.
