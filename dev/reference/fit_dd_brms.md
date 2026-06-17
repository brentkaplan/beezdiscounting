# Fit a Bayesian Mixed-Effects Discounting Model via brms

Fits the four discounting equations of
[`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)
(`"mazur"`, `"exponential"`, `"green-myerson"`, `"rachlin"`) as Bayesian
nonlinear mixed-effects models, with the discount rate estimated on the
natural-log scale (`logk`, subject random intercept; `k ~ 1`, the v1
scope) and the shape exponent (`logs`) population-level for the
two-parameter equations.

## Usage

``` r
fit_dd_brms(
  data,
  y_var = "y",
  x_var = "x",
  id_var = "id",
  equation = c("mazur", "exponential", "green-myerson", "rachlin"),
  family = c("beta", "gaussian"),
  boundary = c("squeeze", "zoib", "error"),
  random_effects = k ~ 1,
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
  ll = NULL,
  response_scale = c("proportion", "percent", "amount"),
  prior = NULL,
  autoscale_priors = TRUE,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  thin = 1,
  cores = getOption("mc.cores", 1L),
  seed = NA,
  backend = getOption("brms.backend", "rstan"),
  control = list(adapt_delta = 0.95),
  init = c("prior_center", "tmb", "random", "0"),
  sample_prior = "no",
  loo = TRUE,
  file = NULL,
  file_refit = getOption("brms.file_refit", "on_change"),
  verbose = 1,
  ...
)
```

## Arguments

- data:

  Long-format data frame (one row per subject-delay).

- y_var, x_var, id_var:

  Column names (canonical defaults), as in
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).

- equation:

  Discounting equation.

- family:

  `"beta"` or `"gaussian"`; `"sltb"` errors with guidance.

- boundary:

  Boundary handling for the beta family (see Details).

- random_effects:

  v1 supports `k ~ 1` only; `k + phi ~ 1` (supported by the TMB tier)
  errors with the planned brms route.

- factors, factor_interaction, continuous_covariates:

  Fixed-effect design on `logk`, as in
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).

- ll, response_scale:

  Response coercion, as in
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).

- prior:

  Optional `brmsprior`; user rows override the defaults.

- autoscale_priors:

  Anchor the `logk` prior to the median delay (see
  [`default_dd_priors()`](https://brentkaplan.github.io/beezdiscounting/reference/default_dd_priors.md)).

- chains, iter, warmup, thin, cores, seed, backend, control,
  sample_prior:

  MCMC settings passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- init:

  `"prior_center"` (default), `"tmb"` (a quiet
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)
  pre-fit supplies the centers, with prior_center fallback on failure;
  the beta family maps to the TMB sltb pre-fit), `"random"`, or `"0"`;
  or a list/function passed through to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- loo:

  Compute and store
  [`brms::loo()`](https://mc-stan.org/loo/reference/loo.html) at fit
  time.

- file, file_refit:

  Passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) for
  fit caching.

- verbose:

  0 (silent), 1 (messages), 2 (full Stan output).

- ...:

  Passed through to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

An object of class `beezdiscounting_brms`: `model` (posterior
medians/SDs on the estimation scale under TMB names `beta_k`/`log_s`,
plus `variance_components`), `brmsfit`, `subject_pars` (`id`, `k`,
`k_lower`, `k_upper`), `converged` (Rhat \< 1.01, no divergences, bulk
ESS \>= 400), `mcmc_info`, `loo`, `data`, `param_info`, `priors`,
`autoscale_info`.

## Details

`family = "beta"` (default) uses `Beta(link = "identity")` with the mean
squished into `(1e-6, 1 - 1e-6)` – the closest brms analog of the TMB
SLT-beta (`family = "sltb"` has no brms equivalent and errors with this
pointer). Boundary observations (`y` exactly 0 or 1) are handled per
`boundary`: `"squeeze"` (default) applies the Smithson-Verkuilen
transform `y* = (y (N - 1) + 0.5) / N` to all responses (message reports
the boundary count); `"zoib"` switches to `zero_one_inflated_beta`
(statistically more honest but changes the estimand – k then describes
interior responses only); `"error"` refuses to fit.
`family = "gaussian"` matches `fit_dd_tmb(family = "gaussian")` wherever
the TMB template's mu clamp into `[1e-6, 1 - 1e-6]` does not bind
(everywhere except extreme decay underflow).

## See also

[`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md);
[`default_dd_priors()`](https://brentkaplan.github.io/beezdiscounting/reference/default_dd_priors.md).
