# Fit a Bayesian Structural Choice Discounting Model via brms

The TMB structural likelihood of
[`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md):
`logit P(LL) = [b0] + gamma * ((ll/ss) * D(k, delay) - 1)` with choice
sensitivity `gamma = exp(loggamma)` and the discount rate
`k = exp(logk)` carrying the subject random intercept. With
`bernoulli("logit")` the brms nonlinear formula IS the logit, so the
likelihood matches TMB exactly.

## Usage

``` r
fit_dd_choice_brms(
  data,
  mode = c("structural", "descriptive"),
  id_var = "id",
  ss_var = "ss_amount",
  ll_var = "ll_amount",
  delay_var = "delay",
  choice_var = "choice",
  equation = c("mazur", "exponential"),
  intercept = FALSE,
  factors = NULL,
  factor_interaction = FALSE,
  continuous_covariates = NULL,
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

  Long-format trial-level data (one row per choice).

- mode:

  `"structural"` (only mode in v1).

- id_var, ss_var, ll_var, delay_var, choice_var:

  Column names, as in
  [`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md).

- equation:

  `"mazur"` or `"exponential"`.

- intercept:

  Include the logit-scale bias term `b0`.

- factors, factor_interaction, continuous_covariates:

  Between-subject fixed-effect design on `log k` (same semantics as
  [`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md));
  `gamma` and `b0` stay population-level.

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

An object of class `beezdiscounting_choice_brms`. Coefficients are
posterior medians on the estimation scale under the TMB names (one
`beta_k` per design column, `log_gamma`, and `beta0` when
`intercept = TRUE`).

## Details

v1 implements `mode = "structural"` only: the descriptive (Young 2018)
model is a plain logistic GLMM expressible directly with
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## See also

[`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md);
[`default_dd_choice_priors()`](https://brentkaplan.github.io/beezdiscounting/reference/default_dd_choice_priors.md);
[`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)
and
[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)
for draws-based marginal means and contrasts of `k`.
