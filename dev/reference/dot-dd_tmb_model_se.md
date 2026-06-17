# Named standard-error vector aligned to the coefficient vector

Pulls fixed-effect SEs from the sdreport. The optimizer parameter vector
(`beta_k`, `log_sigma_u`, `log_aux`) is what `sdreport$par.fixed` /
`sdreport$cov.fixed` cover; entries are renamed to match
`names(object$model$coefficients)` (i.e. `log_aux` already displayed as
`log_phi` / `log_sigma_e` in the fit object). Returns `NA` SEs when the
sdreport is unavailable.

## Usage

``` r
.dd_tmb_model_se(object)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

## Value

Named numeric vector parallel to `object$model$coefficients`.
