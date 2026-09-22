# Predict from a structural choice discounting model

Predict from a structural choice discounting model

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
predict(
  object,
  newdata = NULL,
  type = c("prob", "parameters"),
  level = "subject",
  ...
)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

- newdata:

  Optional data frame. `NULL` uses the fitting data. A supplied
  `newdata` must use the package's canonical column names (`id`,
  `ss_amount`, `ll_amount`, `delay`, plus any factor/covariate columns).

- type:

  `"prob"` (default; fitted choice probabilities P(LL) on `(0,1)`) or
  `"parameters"` (the per-subject parameter tibble).

- level:

  For `type = "prob"`: `"subject"` (default; conditions on each
  subject's estimated random intercept, requires the id column) or
  `"population"` (random effects set to zero – the typical subject: in
  structural mode `log k` at its population value, in descriptive mode
  the random slopes at zero; no id column needed). This is not the
  population-averaged choice probability, which integrates over the
  random effects and is flatter.

- ...:

  Unused.

## Value

- `type = "parameters"`: the subject-parameter tibble (one row per
  subject).

- `type = "prob"`: `newdata` as a tibble plus a `.prob` column of fitted
  choice probabilities.
