# Simulate IP-family mixed-effects discounting data

Generates long-format indifference-point data from the mixed-effects
discounting model used by
[`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).
Each subject `i` has a random discount rate
`log k_i = log_k_pop + delta_k[condition] + u_i`,
`u_i ~ N(0, sigma_u^2)`. The mean indifference proportion at delay `x`
is the discounting function `mu` (Mazur hyperbola, exponential,
Green-Myerson, or Rachlin), and observed `y` is drawn from the
scale-location-truncated beta (`family = "sltb"`) via the inverse-CDF on
the truncated beta, or from a clamped Gaussian (`family = "gaussian"`).

## Usage

``` r
simulate_dd_ip(
  n_subjects = 60,
  delays = c(7, 30, 180, 365, 730, 1460, 2920),
  log_k_pop = log(0.01),
  sigma_u = 0.6,
  phi = 10,
  sigma_e = 0.1,
  family = c("sltb", "gaussian"),
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  s = 1,
  n_conditions = 1,
  delta_k = NULL,
  seed = NULL,
  sigma_phi = 0,
  rho_kphi = 0,
  attach_truth = FALSE,
  sigma_s = 0,
  rho_ks = 0
)
```

## Arguments

- n_subjects:

  Integer; number of subjects.

- delays:

  Numeric vector of delays (days) each subject is observed at.

- log_k_pop:

  Numeric; population intercept on `log k`.

- sigma_u:

  Numeric; SD of the subject random intercept on `log k`.

- phi:

  Numeric; SLT-beta precision (`family = "sltb"`).

- sigma_e:

  Numeric; residual SD on `y` (`family = "gaussian"`).

- family:

  One of `"sltb"` (default) or `"gaussian"`.

- equation:

  One of `"mazur"` (default), `"exponential"`, `"green-myerson"`, or
  `"rachlin"`. The two 2-parameter forms use the nonlinearity exponent
  `s` and reduce to `"mazur"` at `s = 1`.

- s:

  Numeric nonlinearity exponent for the 2-parameter equations
  (Green-Myerson / Rachlin); ignored by `"mazur"` / `"exponential"`.

- n_conditions:

  Integer; number of between-subject condition levels. When `> 1`, a
  `condition` factor is added and subjects are split across levels.

- delta_k:

  Numeric vector of length `n_conditions`; per-condition shift on
  `log k` (the first element is typically `0` for the reference level).
  Required (non-`NULL`) when `n_conditions > 1`.

- seed:

  Optional integer seed.

- sigma_phi:

  Numeric; SD of the subject random intercept on `log phi` (SLT-beta
  precision). Default `0` (no subject-random phi, `family = "sltb"`
  only). When `> 0`, the `(log k, log phi)` pair is drawn jointly from
  `N(0, Sigma)` where `Sigma` is parameterized by `sigma_u`,
  `sigma_phi`, and `rho_kphi`. Requires `family = "sltb"`.

- rho_kphi:

  Numeric in `[-1, 1]`; correlation between subject random intercepts on
  `log k` and `log phi`. Ignored when `sigma_phi = 0`.

- attach_truth:

  Logical; when `TRUE`, a `phi` column is always appended (equal to the
  constant `phi` when `sigma_phi = 0`, or subject-specific when
  `sigma_phi > 0`); when `sigma_s > 0`, an `s` column is also appended
  with the subject-specific curvature exponent. Default `FALSE`.

- sigma_s:

  Numeric; SD of the subject random intercept on `log s` (GM/Rachlin
  curvature exponent). Default `0` (no subject-random s). When `> 0`,
  the `(log k, log s)` pair is drawn jointly from `N(0, Sigma)` where
  `Sigma` is parameterized by `sigma_u`, `sigma_s`, and `rho_ks`.
  Requires `equation = "green-myerson"` or `"rachlin"`. Mutually
  exclusive with `sigma_phi > 0` (q = 2 only).

- rho_ks:

  Numeric in `[-1, 1]`; correlation between subject random intercepts on
  `log k` and `log s`. Ignored when `sigma_s = 0`.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with
columns `id` (factor), `condition` (factor; only when
`n_conditions > 1`), `x` (delay), and `y` (indifference proportion in
`[0, 1]`). When `attach_truth = TRUE`, an additional `phi` column is
appended with the subject-specific SLT-beta precision (equal to the
constant `phi` when `sigma_phi = 0`). When `attach_truth = TRUE` and
`sigma_s > 0`, an additional `s` column is appended with the
subject-specific curvature exponent.

## Details

The SLT draw uses the same constants as the C++ template and the
verified reference density: `s_slt = 1.0000001`, `l = 1e-8`, with
`a = mu * phi`, `b = (1 - mu) * phi`, and
`y = (qbeta(U, a, b) - l) * s_slt` where
`U ~ Uniform(pbeta(l, a, b), pbeta(1/s_slt + l, a, b))`.

## Examples

``` r
# A small SLT-beta mixed-effects discounting dataset (id, x, y)
sim <- simulate_dd_ip(n_subjects = 8, seed = 1)
head(sim)
#> # A tibble: 6 × 3
#>   id        x      y
#>   <fct> <dbl>  <dbl>
#> 1 1         7 0.995 
#> 2 1        30 0.991 
#> 3 1       180 0.395 
#> 4 1       365 0.388 
#> 5 1       730 0.361 
#> 6 1      1460 0.0210

# Two between-subject conditions with a log-k shift on the second level
sim2 <- simulate_dd_ip(
  n_subjects = 8, n_conditions = 2, delta_k = c(0, log(3)), seed = 2
)
```
