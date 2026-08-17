# Tidy a beezdiscounting_choice model into a coefficient tibble

Returns fixed-effect rows (the log-k coefficients), the shape rows (the
choice-sensitivity `gamma` and, when present, the choice-bias `beta0`),
and the variance-component row, following the broom coefficient-table
contract.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
tidy(
  x,
  effects = c("fixed", "ran_pars"),
  report_space = c("natural", "log10", "internal", "log"),
  ...
)
```

## Arguments

- x:

  A `beezdiscounting_choice` object.

- effects:

  Character vector: `"fixed"` (log-k fixed-effect rows + the shape rows
  gamma/beta0), `"ran_pars"` (the RE SD), or both (default).

- report_space:

  `"natural"`, `"log10"`, `"internal"`, or `"log"` - reporting scale for
  the fixed-effect `estimate`/`std.error`. Default is `"natural"`.

- ...:

  Unused.

## Value

A tibble with exactly 8 columns in this order: `term`, `estimate`,
`std.error`, `statistic`, `p.value`, `component`, `estimate_scale`,
`term_display`. Fixed-effect rows carry `component == "fixed"`;
gamma/beta0 carry `component == "shape"`; variance rows carry
`component == "variance"`.

## Details

`estimate` and `std.error` are reported on the `report_space` scale for
the fixed-effect (`beta_k`) rows and for `gamma` (which is transformed
EXPLICITLY since the param-space transformer keys only on k/s/phi).
`beta0` is on the identity (logit-intercept) scale and is NEVER
transformed across report spaces. `statistic` and `p.value` are always
computed on the estimation (internal) scale - Wald statistics are not
recomputed after back-transforming (broom convention).
