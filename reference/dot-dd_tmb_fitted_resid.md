# Shared fitted-values / residuals back end

Calls [`predict()`](https://rdrr.io/r/stats/predict.html) at the
requested level and returns `.fitted`, `.resid`, and the data frame used
(either `object$data` or `newdata`).

## Usage

``` r
.dd_tmb_fitted_resid(
  object,
  newdata = NULL,
  level = c("subject", "population")
)
```

## Arguments

- object:

  A `beezdiscounting_tmb` fit.

- newdata:

  Optional data frame; `NULL` -\> training data.

- level:

  `"subject"` (default) or `"population"`.

## Value

List with `.fitted`, `.resid`, and `data`.
