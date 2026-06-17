# Build display term names from a beezdiscounting_tmb coefficient vector

Maps raw optimizer names (`beta_k`, `log_sigma_u`, `log_phi` /
`log_sigma_e`) to readable display names. `beta_k` columns become
`k:<design colname>` using `formula_details$X` colnames; every other
coefficient keeps its raw name. Shared by
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint()`](https://rdrr.io/r/stats/confint.html).

## Usage

``` r
.dd_tmb_build_term_names(object, nms = NULL)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- nms:

  Character vector of raw parameter names (default
  `names(object$model$coefficients)`).

## Value

List with `term` (display names), `k_idx` (beta_k positions),
`other_idx` (non-beta positions).
