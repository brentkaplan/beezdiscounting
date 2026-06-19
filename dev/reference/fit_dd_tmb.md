# Fit an indifference-point mixed-effects discounting model via TMB

Fits a 1-parameter discounting model (Mazur hyperbolic or exponential)
with a random intercept on `log k`, between-subject fixed effects, and
either an SLT-beta or Gaussian observation family, using Template Model
Builder for exact AD + Laplace approximation.

## Usage

``` r
fit_dd_tmb(
  data,
  y_var = "y",
  x_var = "x",
  id_var = "id",
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  family = c("sltb", "gaussian"),
  random_effects = k ~ 1,
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  ll = NULL,
  response_scale = c("proportion", "percent", "amount"),
  start_values = NULL,
  tmb_control = list(iter_max = 1000, eval_max = 2000),
  multi_start = TRUE,
  verbose = 1,
  covariance_structure = c("pdSymm", "pdDiag"),
  ...
)
```

## Arguments

- data:

  Long data frame with subject id, delay, and indifference proportion
  columns.

- y_var, x_var, id_var:

  Column names (defaults `"y"`, `"x"`, `"id"`).

- equation:

  One of `"mazur"`, `"exponential"`, `"green-myerson"`, or `"rachlin"`.
  The two 2-parameter (hyperboloid) forms add a single population
  nonlinearity exponent `s` (estimated on the log scale) and reduce to
  `"mazur"` at `s = 1`.

- family:

  Observation family: `"sltb"` (default) or `"gaussian"`.

- random_effects:

  RE formula: `k ~ 1` (single random intercept on `log k`),
  `k + phi ~ 1` (a joint 2-D random intercept on `(log k, log phi)`,
  SLT-beta only), or `k + s ~ 1` (a joint 2-D random intercept on
  `(log k, log s)`, Green-Myerson and Rachlin only).

- factors:

  Character vector of between-subject factor names.

- factor_interaction:

  Logical; include a pairwise factor interaction.

- continuous_covariates:

  Character vector of covariate names.

- ll:

  Optional larger-later reward for `amount`-scale coercion.

- response_scale:

  One of `"proportion"`, `"percent"`, `"amount"`.

- start_values:

  Optional named list overriding defaults.

- tmb_control:

  Optimizer control list.

- multi_start:

  Logical; if `TRUE` (default), run the 3-set guarded multi-start.

- verbose:

  Integer verbosity (0 silent, 1 progress, 2 debug).

- covariance_structure:

  Covariance for a 2-D random effect (`k + phi ~ 1` or `k + s ~ 1`):
  `"pdSymm"` (default; correlated random intercepts) or `"pdDiag"`
  (independent, correlation fixed at 0). Ignored for `k ~ 1`.

- ...:

  Reserved.

## Value

An object of class `beezdiscounting_tmb` with components:

- call:

  The matched call.

- opt:

  Normalized optimizer result (`par`, `objective`, `convergence`,
  `message`).

- model:

  List of `coefficients`, `se`, and `variance_components`.

- sdr:

  TMB `sdreport` object (or `NULL` if SE computation failed).

- hessian_pd:

  Logical positive-definiteness of the Hessian.

- param_info:

  Model metadata (equation, family, dimensions, factor spec, parsed
  random effects).

- formula_details:

  Fixed-effect design (`X`, `rhs`, `contrasts`).

- subject_pars:

  Data frame of subject-level parameters. For a 1-RE fit (`k ~ 1`) the
  columns are `id, u_i, k`; for a phi-target 2-RE fit (`k + phi ~ 1`)
  they are `id, re_k, re_phi, k, phi`; for an s-target 2-RE fit
  (`k + s ~ 1`, GM/Rachlin) they are `id, re_k, re_s, k, s` where `s` is
  soft-clamped toward `(0.05, 20)`.

- loglik, AIC, BIC:

  Fit statistics.

- converged, se_available:

  Convergence / SE-availability flags.

- opt_warnings:

  Character vector of optimizer warnings.

- data:

  The single filtered model frame (id/x/y + retained design columns),
  row-aligned with the design matrix.

- data_all:

  The validated frame before complete-casing.

- coercion_info:

  Scale-coercion/clamping audit list.

## References

Young, M. E. (2017). Discounting: A practical guide to multilevel
analysis of indifference data. *Journal of the Experimental Analysis of
Behavior, 108*(1), 97-112.
[doi:10.1002/jeab.265](https://doi.org/10.1002/jeab.265)

Kim, M., Koffarnus, M. N., & Franck, C. T. (2024). Thinking inside the
bounds: Improved error distributions for indifference point data
analysis and simulation via beta regression using common discounting
functions. *arXiv* preprint arXiv:2404.18000.

Kim, M., Kaplan, B. A., Koffarnus, M. N., & Franck, C. T. (2025).
Scale-location-truncated beta regression: Expanding beta regression to
accommodate 0 and 1. *arXiv* preprint arXiv:2509.13167.

## Examples

``` r
# \donttest{
# Small two-subject long-format indifference-point data frame.
dd <- data.frame(
  id = rep(c("s1", "s2"), each = 5),
  x  = rep(c(7, 30, 180, 365, 730), times = 2),
  y  = c(0.95, 0.80, 0.45, 0.30, 0.15,
         0.90, 0.70, 0.40, 0.25, 0.10)
)
fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb", verbose = 0)
exp(fit$model$coefficients[["beta_k"]])  # population k
#> [1] 0.008690582
# }
```
