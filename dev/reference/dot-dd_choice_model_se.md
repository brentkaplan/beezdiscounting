# Named standard-error vector aligned to the coefficient vector

Returns `object$model$se` aligned to `names(object$model$coefficients)`.
Returns `NA` SEs when `se_available` is `FALSE` (non-PD Hessian or
failed sdreport), so confint()/tidy()/summary() never present unreliable
Wald intervals/p-values.

## Usage

``` r
.dd_choice_model_se(object)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

## Value

Named numeric vector parallel to `object$model$coefficients`.
