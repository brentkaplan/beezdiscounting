# Find the smallest sample size reaching a target power (discounting)

Bisection search over `n_subjects` for the smallest total N whose Monte
Carlo power estimate from
[`power_discounting()`](https://brentkaplan.github.io/beezdiscounting/reference/power_discounting.md)
reaches `target_power`. The search accounts for Monte Carlo noise: at
each evaluated N, replicates are added in batches (up to `n_sim_max`)
until the Wilson interval for power lies wholly above or below the
target; if it still straddles the target at the cap, the decision falls
back to the point estimate and the result is flagged `uncertain`. The
selected N and its lower neighbor are then re-evaluated with fresh
replicates before minimality is claimed.

The returned `n` is an *estimated minimum under Monte Carlo uncertainty*
rather than an exact bound. For grant-quality reporting, rerun
[`power_discounting()`](https://brentkaplan.github.io/beezdiscounting/reference/power_discounting.md)
at the returned `n` with a large `n_sim` (2000+) and report that
estimate with its Monte Carlo confidence interval.

**Monotonicity assumption.** Bisection presumes that power is
non-decreasing in `n_subjects`. Because every evaluated N is judged from
its own independent replicates (and a convergence-conditioned
denominator), a Monte Carlo fluctuation can make a lower N read "below"
when its true power is above the target, so the search may step past a
lower crossing that it never revisits. Evaluated N that contradict the
assumption (a lower N reading "above" the selected N) demote the status
to `"uncertain"`; N that were never evaluated cannot be checked. Widen
`n_sim`/`n_sim_max` when the reported `n` matters. When the target is
already met at `n_range[1]`, that bound is likewise re-evaluated with
fresh replicates before `"at_lower_bound"` is reported; if the second
look does not clear the target the bound is treated as below and the
bisection proceeds upward.

## Usage

``` r
find_n_discounting(
  target_power = 0.8,
  effect = list(delta_k = NULL),
  design = list(),
  n_range = c(6, 200),
  n_sim = 200,
  n_sim_max = 4 * n_sim,
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

- target_power:

  Target power in (0, 1).

- effect:

  Named list supplying `delta_k`: the true condition-2 shift on
  natural-log k. `0` is allowed (useful for Type I error checks). E.g.
  `delta_k = log(2)` means condition 2's k is twice condition 1's.

- design:

  Named list of data-generating settings, merged over the simulator
  defaults: `delays` (vector), `log_k_pop`, `sigma_u` (subject SD on log
  k), `phi` (SLT-beta precision; used when `family = "sltb"`), and
  `sigma_e` (residual SD; used when `family = "gaussian"`).

- n_range:

  Integer bracket `c(lower, upper)` of total N to search
  (`2 <= lower < upper`). The search errors (rather than extrapolating)
  if the target is not reached at `upper`.

- n_sim:

  Replicates per evaluation batch. Smaller than the
  [`power_discounting()`](https://brentkaplan.github.io/beezdiscounting/reference/power_discounting.md)
  default because several N values are evaluated; the adaptive rule adds
  batches where the verdict is close.

- n_sim_max:

  Maximum replicates per evaluated N (default `4 * n_sim`).

- alpha:

  Nominal two-sided test level.

- df:

  Degrees of freedom for the Wald test's t reference. `NULL` (default)
  tracks the evaluated sample size as `n - 2`; a numeric value (or `Inf`
  for the asymptotic z-test) is used at every evaluated N.

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
  refit of data with no such variance component. That refit is useful
  for probing robustness rather than for estimating power under those
  random effects (out of scope in this version).

- multi_start:

  Passed to
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).
  Defaults to `FALSE` for speed; non-convergent replicates are excluded
  and surfaced rather than biasing the estimate.

- verbose:

  Logical; report each evaluation.

- ...:

  Additional arguments passed to
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).

## Value

An object of class `beezdiscounting_power_n`: a list with

- n:

  Estimated smallest total `n_subjects` reaching `target_power`; `NA`
  when the confirmation pass contradicted the search
  (`status = "unresolved"`).

- target_power:

  As supplied.

- status:

  `"confirmed"` (selected N re-confirmed above target and N - 1 below),
  `"uncertain"` (a decision relied on a point estimate, N - 1 also
  cleared the target on reconfirmation, or an evaluated lower N read
  above the target, so the returned N may not be minimal),
  `"unresolved"` (the selected N failed its own reconfirmation; `n` is
  `NA`), or `"at_lower_bound"` (the target was already met at
  `n_range[1]` on two independent looks; smaller N was not explored;
  widen `n_range` downward if that matters). These labels describe a
  heuristic Monte Carlo decision rule (repeated looks at ordinary Wilson
  intervals across several N) rather than a formal sequential error
  guarantee.

- uncertain:

  Logical; `TRUE` when any search decision was made on a point estimate
  rather than a conclusive Wilson interval, or the status is not
  `"confirmed"`/`"at_lower_bound"`.

- evaluations:

  Tibble of every evaluation: `n_subjects`, `n_sim_total`, `n_used`,
  `usable_fraction`, `power`, `ci_lower`, `ci_upper`, `decision`. A
  warning fires when any evaluation had fewer than 95% usable fits.

- alpha, df, effect, design, n_range, n_sim, n_sim_max, seed, settings,
  call:

  Echoed inputs and effective settings.

## See also

[`power_discounting()`](https://brentkaplan.github.io/beezdiscounting/reference/power_discounting.md)
for the Monte Carlo engine.

Other power-analysis:
[`power_discounting()`](https://brentkaplan.github.io/beezdiscounting/reference/power_discounting.md)

## Examples

``` r
# \donttest{
# Small search for demonstration (use larger n_sim for real planning)
res <- find_n_discounting(
  target_power = 0.8,
  effect = list(delta_k = log(8)),
  n_range = c(6, 20), n_sim = 4, n_sim_max = 4, seed = 1, verbose = FALSE
)
print(res)
#> Sample-size search (Monte Carlo power)
#>   Target power 0.80 for delta_k = 2.079 at alpha = 0.05
#>   Estimated minimum n_subjects = 8 (status: uncertain)
#>   This is an estimated minimum under Monte Carlo uncertainty;
#>   rerun the power function at this N with a large n_sim to report it.
#> 
#>   Evaluations:
#>  n_subjects n_sim_total n_used usable_fraction power  ci_lower  ci_upper
#>          20           4      4               1  1.00 0.5101092 1.0000000
#>           6           4      4               1  0.75 0.3006418 0.9544127
#>          13           4      4               1  1.00 0.5101092 1.0000000
#>           9           4      4               1  1.00 0.5101092 1.0000000
#>           7           4      4               1  0.75 0.3006418 0.9544127
#>           8           4      4               1  1.00 0.5101092 1.0000000
#>           8           4      4               1  1.00 0.5101092 1.0000000
#>           7           4      4               1  1.00 0.5101092 1.0000000
#>         decision
#>  ambiguous_above
#>  ambiguous_below
#>  ambiguous_above
#>  ambiguous_above
#>  ambiguous_below
#>  ambiguous_above
#>  ambiguous_above
#>  ambiguous_above
# }
```
