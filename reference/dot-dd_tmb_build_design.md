# Build the fixed-effect design matrix for log k

Constructs the `model.matrix` for the `log k = Xbeta` linear predictor
from between-subject factors, an optional pairwise interaction, and
continuous covariates. With no inputs the design is intercept-only.

## Usage

``` r
.dd_tmb_build_design(
  data,
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL
)
```

## Arguments

- data:

  Cleaned long data frame (canonical `id`/`x`/`y` plus any factor or
  covariate columns).

- factors:

  Character vector of factor column names, or `NULL`.

- factor_interaction:

  Logical; if `TRUE` and \>= 2 factors, include their interaction (uses
  `*`); otherwise main effects (`+`).

- continuous_covariates:

  Character vector of covariate names, or `NULL`.

## Value

A list with `X` (the model matrix for log-k FE), `rhs` (the one-sided
RHS **formula** from
[`build_fixed_rhs()`](https://brentkaplan.github.io/beezdemand/reference/build_fixed_rhs.html)),
and `contrasts` (`attr(X, "contrasts")`, the per-factor contrast
scheme). The `rhs` and `contrasts` are stored so
`.dd_build_emm_ref_grid()` can rebuild a column-aligned design via the
same route.
