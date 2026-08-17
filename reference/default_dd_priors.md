# Default priors for Bayesian (brms) delay-discounting models

Returns the default prior table used by
[`fit_dd_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_brms.md),
so the defaults can be inspected, modified row-wise, and passed back via
the fitter's `prior` argument. The `logk` location is the principled fix
for k's delay-unit dependence: with `autoscale = TRUE` (the default
whenever `data` is supplied) it centers `k * median(delay) = 1` – the
delay at which the Mazur curve crosses 0.5 – via
`normal(-log(median(x)), 2.5)`; the static fallback is
`normal(-4.5, 2.5)`. The anchors used are attached as
`attr(, "autoscale_info")`; numeric values are formatted with
`format(x, digits = 6, scientific = FALSE)`.

## Usage

``` r
default_dd_priors(
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  family = c("beta", "gaussian"),
  data = NULL,
  y_var = "y",
  x_var = "x",
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  autoscale = !is.null(data),
  random_effects = k ~ 1,
  covariance_structure = c("pdSymm", "pdDiag")
)
```

## Arguments

- equation:

  Discounting equation (TMB-tier vocabulary).

- family:

  `"beta"` or `"gaussian"`.

- data:

  Optional data frame used for autoscaling.

- y_var, x_var:

  Column names in `data` (canonical defaults).

- factors, factor_interaction, continuous_covariates:

  Fixed-effect design on `logk`, as passed to the fitter. When the
  design has non-intercept coefficients (derived through
  `build_fixed_rhs()`, so single-level dropped factors do not count), a
  fold-change `normal(0, 1)` class-level coefficient prior is added;
  with an intercept-only design it is omitted (it would be unused, and
  brms warns).

- autoscale:

  Logical; defaults to `TRUE` when `data` is supplied.

- random_effects:

  `k ~ 1` (default) or `k + phi ~ 1`. The latter re-keys the beta
  precision: phi becomes a predicted distributional parameter, so the
  scalar `gamma(2, 0.1)` is replaced by a log-scale intercept prior, a
  half-t precision-RE SD, and (for `"pdSymm"`) an LKJ correlation prior.

- covariance_structure:

  `(log k, log phi)` covariance for `k + phi ~ 1`: `"pdSymm"` (default,
  correlated) or `"pdDiag"` (independent).

## Value

A `brmsprior` data frame, with `attr(, "autoscale_info")` when
autoscaling was used.

## Details

Other defaults: `logs ~ normal(0, 0.5)` (two-parameter equations; s is
near 1 a priori), `sd(logk) ~ student_t(3, 0, 1)`, `phi ~ gamma(2, 0.1)`
(Beta precision; mean 20, far from brms's near-improper
`gamma(0.01, 0.01)`), and `sigma ~ student_t(3, 0, 0.25)` for the
Gaussian family (y is a proportion in the unit interval, so sd(y) \<=
0.5).
