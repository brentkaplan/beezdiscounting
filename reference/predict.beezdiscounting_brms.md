# Predict from a beezdiscounting_brms model

`type = "response"` summarizes
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
draws (the indifference-proportion mean) to the posterior median with
equal-tailed intervals; `level = "subject"` conditions on the subject
random effects, `"population"` sets them to zero. `type = "parameters"`
returns the per-subject k summaries.

## Usage

``` r
# S3 method for class 'beezdiscounting_brms'
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

  A `beezdiscounting_brms` object.

- newdata:

  Optional data frame (defaults to the model frame). The Rachlin guard
  columns (`xzero`/`xsafe`) are derived automatically.

- type:

  `"response"` or `"parameters"`.

- level:

  `"subject"` (default) or `"population"`.

- probs:

  Interval probabilities.

- draws:

  If `TRUE`, return the epred draws matrix.

- ...:

  Unused.

## Value

A tibble (or draws matrix).
