# Default priors for the Bayesian (brms) choice model

`loggamma ~ normal(1, 1)` (choice sensitivity, typically 1-20) and the
optional logit-scale bias `b0 ~ normal(0, 1.5)`; `logk` defaults as in
[`default_dd_priors()`](https://brentkaplan.github.io/beezdiscounting/reference/default_dd_priors.md)
(autoscaled to the median delay when `data` is supplied; `delay_var`
names the delay column).

## Usage

``` r
default_dd_choice_priors(
  equation = c("mazur", "exponential"),
  intercept = FALSE,
  data = NULL,
  delay_var = "delay",
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  autoscale = !is.null(data)
)
```

## Arguments

- equation:

  `"mazur"` or `"exponential"`.

- intercept:

  Include the `b0` prior.

- data:

  Optional data frame used for `logk` autoscaling.

- delay_var:

  Delay column name in `data`.

- factors, factor_interaction, continuous_covariates:

  Fixed-effect design on `logk`, as passed to
  [`fit_dd_choice_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice_brms.md).
  When the design has non-intercept coefficients (derived through
  `build_fixed_rhs()`), a fold-change `normal(0, 1)` class-level
  coefficient prior is added; with an intercept-only design it is
  omitted (it would be unused, and brms warns).

- autoscale:

  Logical; defaults to `TRUE` when `data` is supplied.

## Value

A `brmsprior` data frame.
