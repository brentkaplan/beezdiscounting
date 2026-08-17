# Tidy a beezdiscounting_tmb model into a coefficient tibble

Returns fixed-effect rows (the log-k coefficients) and/or
variance-component rows, following the broom coefficient-table contract.

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
tidy(
  x,
  effects = c("fixed", "ran_pars"),
  report_space = c("natural", "log10", "internal", "log"),
  ...
)
```

## Arguments

- x:

  A `beezdiscounting_tmb` object.

- effects:

  Character vector: `"fixed"` (log-k fixed-effect rows), `"ran_pars"`
  (the RE SD and the auxiliary precision/scale parameter), or both
  (default).

- report_space:

  `"natural"`, `"log10"`, `"internal"`, or `"log"` - reporting scale for
  fixed-effect `estimate`/`std.error` (`"internal"` and `"log"` coincide
  for the log-k `beta_k` rows). Default is `"natural"`.

- ...:

  Unused.

## Value

A tibble with exactly 8 columns in this order: `term`, `estimate`,
`std.error`, `statistic`, `p.value`, `component`, `estimate_scale`,
`term_display`. Fixed-effect rows carry `component == "fixed"`; the
shape-parameter row (2-parameter equations only) carries
`component == "shape"`; variance rows carry `component == "variance"`.

## Details

`estimate` and `std.error` are reported on the `report_space` scale for
the fixed-effect (`beta_k`) rows. `statistic` and `p.value` are always
computed on the estimation (log-k) scale - Wald statistics are not
recomputed after back-transforming (broom convention; see the
[`summary()`](https://rdrr.io/r/base/summary.html) note for details).
Variance-component rows carry `NA` for `statistic` and `p.value` and are
not affected by `report_space`.
