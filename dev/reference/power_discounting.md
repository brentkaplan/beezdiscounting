# Monte Carlo power analysis for between-subject discounting designs

Estimates statistical power to detect a between-subject difference in
the discount rate (a shift `delta_k` on natural-log k) by simulation:
each replicate (1) simulates a two-condition between-subject
indifference-point dataset with
[`simulate_dd_ip()`](https://brentkaplan.github.io/beezdiscounting/reference/simulate_dd_ip.md)
under assumed population parameters plus the effect, (2) refits it with
[`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md),
and (3) tests the log-k condition contrast with a Wald test at level
`alpha`, referred to a t distribution with `df` degrees of freedom (see
the `df` argument). Power is the proportion of *usable* fits (converged,
positive-definite Hessian, finite standard error) that reject.

Because the power estimate is a proportion from finitely many
replicates, it is reported with a Wilson score confidence interval
(`power_mc_ci`). Both a p-value verdict (`p < alpha`) and a
confidence-interval verdict (Wald CI excludes 0) are recorded per
replicate; they use the same standard error and reference distribution,
so they coincide by construction, and both rates are returned.

This mirrors `beezdemand::power_demand()`. The design asymmetry is
intentional: this package's simulator models a *between-subject*
condition (subjects are split across conditions round-robin), while
beezdemand's models a within-subject condition.

## Usage

``` r
power_discounting(
  n_subjects,
  effect = list(delta_k = NULL),
  design = list(),
  n_sim = 500,
  alpha = 0.05,
  df = NULL,
  seed = NULL,
  equation = c("mazur", "exponential"),
  family = c("sltb", "gaussian"),
  random_effects = k ~ 1,
  multi_start = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- n_subjects:

  Total number of simulated subjects per replicate, split across the two
  conditions round-robin (even numbers give equal groups).

- effect:

  Named list supplying `delta_k`: the true condition-2 shift on
  natural-log k. `0` is allowed (useful for Type I error checks). E.g.
  `delta_k = log(2)` means condition 2's k is twice condition 1's.

- design:

  Named list of data-generating settings, merged over the simulator
  defaults: `delays` (vector), `log_k_pop`, `sigma_u` (subject SD on log
  k), `phi` (SLT-beta precision; used when `family = "sltb"`), and
  `sigma_e` (residual SD; used when `family = "gaussian"`).

- n_sim:

  Number of Monte Carlo replicates. 500 (default) is suitable for
  interactive exploration; use 2000+ for grant-quality precision (see
  [`vignette("power-analysis")`](https://brentkaplan.github.io/beezdiscounting/articles/power-analysis.md)).

- alpha:

  Nominal two-sided test level.

- df:

  Degrees of freedom for the Wald test's t reference distribution.
  `NULL` (default) uses `n_subjects - 2`, the two-sample df of the
  between-subject design. This is an *empirically calibrated*
  small-sample correction, not a model-derived df (the TMB fit has no
  exact t sampling theory); it passes the package's Type I calibration
  battery, while the asymptotic z-test (`df = Inf`) is anticonservative
  at study-relevant sample sizes.

- seed:

  Optional integer seed; identical seeds give identical results. The
  caller's RNG state is restored on exit.

- equation:

  Discounting function used for BOTH simulation and refitting (they are
  always matched): `"mazur"` (default) or `"exponential"`. The
  two-parameter Green-Myerson and Rachlin forms are out of scope in this
  version.

- family:

  Observation family used for both simulation and refitting: `"sltb"`
  (default) or `"gaussian"`.

- random_effects:

  Random-effects formula passed to
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).
  The default `k ~ 1` matches the simulator's data-generating process (a
  single subject random intercept on log k). The simulator ALWAYS
  generates only that intercept: any richer formula (e.g. `k + phi ~ 1`,
  `k + s ~ 1`) is accepted but produces a deliberately over-specified
  refit of data with no such variance component – useful for probing
  robustness, not for estimating power under those random effects (out
  of scope in this version).

- multi_start:

  Passed to
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).
  Defaults to `FALSE` for speed; non-convergent replicates are excluded
  and surfaced rather than biasing the estimate.

- verbose:

  Logical; show a progress bar.

- ...:

  Additional arguments passed to
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)
  (e.g. `tmb_control`).

## Value

An object of class `beezdiscounting_power`: a list with

- power:

  Estimated power: proportion of usable replicates whose Wald CI
  excludes 0 (equal to `hit_rate_ci`). `NA` if no replicate was usable.

- power_mc_ci:

  Wilson 95% confidence interval on `power`, reflecting Monte Carlo
  uncertainty from `n_used` replicates.

- hit_rate_p:

  Proportion of usable replicates with `p < alpha`.

- hit_rate_ci:

  Proportion of usable replicates whose Wald CI excludes 0 (the same
  decision rule as `hit_rate_p`, since both use the same SE and t
  reference; both reported).

- n_sim:

  Total replicates attempted.

- n_converged:

  Replicates whose fit converged.

- n_hessian_pd:

  Replicates with a positive-definite Hessian.

- n_used:

  Replicates entering the power denominator (converged,
  positive-definite Hessian, finite SE).

- alpha:

  Nominal test level.

- df:

  Degrees of freedom of the t reference distribution actually used
  (`n_subjects - 2` unless overridden).

- effect:

  The validated effect specification (name and delta).

- target_term:

  The tested coefficient (`"k:conditionC2"`).

- design:

  The merged design list actually used.

- n_subjects:

  As supplied.

- replicates:

  Tibble with one row per replicate: `sim`, `status` (`"ok"`,
  `"nonconverged"`, `"hessian_not_pd"`, `"se_unusable"`, `"error"`),
  `converged`, `hessian_pd`, `estimate`, `se`, `statistic`, `p_value`,
  `ci_lower`, `ci_upper`, `hit_p`, `hit_ci`, and `message` (error text,
  if any). Estimates are on the natural-log-k scale of the simulated
  `delta_k`.

- seed:

  As supplied.

- settings:

  List of `equation`, `family`, `multi_start`, and the deparsed
  random-effects formula.

- call:

  The matched call.

## Details

A replicate whose fit fails (non-convergence, non-positive-definite
Hessian, unusable standard error, or an error) is excluded from the
power denominator and reported through the `n_*` counts and
`$replicates$status` – it is never counted as "no effect detected",
which would bias power in an unpredictable direction. A warning is
issued when fewer than 95% of replicates are usable.

The v1 scope is a single fixed-effect `delta_k` under the package's
existing between-subject simulator. Effects on `s` or `phi`, power for
derived measures (ED50, AUC), and arbitrary designs are out of scope;
see
[`vignette("power-analysis")`](https://brentkaplan.github.io/beezdiscounting/articles/power-analysis.md).

## See also

[`find_n_discounting()`](https://brentkaplan.github.io/beezdiscounting/reference/find_n_discounting.md)
to search for the smallest adequate sample size;
[`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)
for the model being refit;
[`simulate_dd_ip()`](https://brentkaplan.github.io/beezdiscounting/reference/simulate_dd_ip.md)
for the data-generating process.

Other power-analysis:
[`find_n_discounting()`](https://brentkaplan.github.io/beezdiscounting/reference/find_n_discounting.md)

## Examples

``` r
# \donttest{
# Quick exploratory run (use n_sim >= 500 for real planning)
res <- power_discounting(
  n_subjects = 8,
  effect = list(delta_k = log(2)),
  n_sim = 4, seed = 1, verbose = FALSE
)
print(res)
#> Monte Carlo power analysis (beezdiscounting)
#>   Target: k:conditionC2 (delta_k = 0.6931), two-sided alpha = 0.05, t reference (df = 6)
#>   n_subjects = 8 (between-subject, 2 conditions), n_sim = 4 (converged 4, usable 4)
#>   Power (CI-exclusion): 0.250 [95% MC CI 0.046, 0.699]
#>   p-value hit rate:     0.250
# }
```
