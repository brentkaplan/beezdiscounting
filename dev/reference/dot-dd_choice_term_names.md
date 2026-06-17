# Build display term names from a beezdiscounting_choice coefficient vector

Maps raw optimizer names to readable display names. `beta_k` columns
become `k:<design colname>` using `formula_details$X` colnames;
`log_gamma` becomes `gamma`; every other coefficient (`log_sigma_u`,
`beta0`) keeps its raw name. Shared by
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`print()`](https://rdrr.io/r/base/print.html).

## Usage

``` r
.dd_choice_term_names(object, nms = NULL)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

- nms:

  Character vector of raw parameter names (default
  `names(object$model$coefficients)`).

## Value

List with `term` (display names), `k_idx` (beta_k positions),
`other_idx` (non-beta positions).
