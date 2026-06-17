# Shared fitted-probabilities / residuals back end

Calls `predict(type = "prob")` at the requested level and returns the
fitted choice probabilities, raw residuals (`choice - p`), and the data
frame used.

## Usage

``` r
.dd_choice_fitted_resid(
  object,
  newdata = NULL,
  level = c("subject", "population")
)
```

## Arguments

- object:

  A `beezdiscounting_choice` fit.

- newdata:

  Optional data frame; `NULL` -\> training data.

- level:

  `"subject"` (default) or `"population"`.

## Value

List with `.fitted`, `.resid`, and `data`.
