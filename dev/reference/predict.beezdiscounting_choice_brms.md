# Predict P(LL) from a beezdiscounting_choice_brms model

Predict P(LL) from a beezdiscounting_choice_brms model

## Usage

``` r
# S3 method for class 'beezdiscounting_choice_brms'
predict(
  object,
  newdata = NULL,
  type = c("response", "parameters"),
  level = "subject",
  probs = c(0.025, 0.975),
  draws = FALSE,
  ...
)
```

## Arguments

- object:

  A `beezdiscounting_choice_brms` object.

- newdata:

  Optional trial-level data (needs `rel` and `delay`, or `ss`/`ll`
  columns from which `rel` is derived).

- type:

  `"response"` (P(choosing the larger-later)) or `"parameters"`
  (per-subject k summaries).

- level:

  `"subject"` or `"population"`.

- probs:

  Interval probabilities.

- draws:

  Return the draws matrix instead.

- ...:

  Unused.

## Value

A tibble (or draws matrix).
