# Per-row response SD on the `[0,1]` scale for standardized residuals

Gaussian: constant `sigma_e = exp(log_sigma_e)`. SLT-beta: the
delta-method SLT SD `s_slt * sqrt(mu * (1 - mu) / (phi + 1))` (the SLT
variance at `s_slt ~ 1`), so residuals near the bounds are
down-weighted.

## Usage

``` r
.dd_tmb_response_sd(object, mu, ids = NULL)
```

## Arguments

- object:

  A `beezdiscounting_tmb` fit.

- mu:

  Numeric vector of fitted mu values (from `.dd_discount_mu`).

- ids:

  Optional per-row subject ids (length `mu`). For a **phi-target** 2-RE
  SLT fit at the subject level, supplying `ids` uses each subject's
  `phi_i` so the SD is subject-conditional. For an **s-target** 2-RE fit
  there is no per-subject `phi`; the SD always uses the population
  precision `exp(log_phi)`. `NULL` (population level, or any 1-RE fit)
  also uses the population precision.

## Value

Numeric vector of per-row response SDs, same length as `mu`.
