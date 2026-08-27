# Fit the linearized Mazur hyperbola (Hinds et al., 2026)

Linearizes `D = 1/(1 + k t)` to `ln(1/D - 1) - ln(t) = ln(k) + e` so
that ln(k) has a closed-form per-subject estimator (the geometric mean
of `(1/D - 1)/t`), a one-way random-effects model on the transformed
scale has closed-form MLEs, and condition means can be compared with an
exact F-test
([`anova.beezdiscounting_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/anova.beezdiscounting_linear.md));
the test is exact when the random-effects variance estimate `g-hat > 0`,
and when `g-hat = 0` the statistic is still reported, with a warning. No
numerical optimization is involved.

## Usage

``` r
fit_dd_linear(
  data,
  y_var = "y",
  x_var = "x",
  id_var = "id",
  factors = NULL,
  response_scale = c("proportion", "percent", "amount"),
  ll = NULL,
  boundary = c("clamp", "drop", "error"),
  eps = NULL,
  conf_level = 0.95
)
```

## Arguments

- data:

  Long data frame with subject id, delay, and indifference point
  columns.

- y_var, x_var, id_var:

  Column names for indifference point, delay, subject.

- factors:

  Optional name of ONE categorical column giving the experimental
  condition. Must be between-subject: every id appears under exactly one
  level; frames with an id under several levels are rejected because the
  exact F-test assumes independent units per condition.

- response_scale, ll:

  As in
  [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md).

- boundary:

  How to treat `y` in {0, 1}: `"clamp"` (default; exact 0 becomes `eps`,
  exact 1 becomes `1 - eps`, nothing else moves), `"drop"` (per-subject
  estimates only unless the design stays balanced), `"error"`.

- eps:

  Where exact 0 / 1 are placed under `boundary = "clamp"` (`eps` and
  `1 - eps`); default `1/(2 * ll)` if `ll` is given, else `0.005`.

- conf_level:

  Confidence level for per-subject ln(k) intervals (a single number
  strictly between 0 and 1); also the default `level` of
  [`confint()`](https://rdrr.io/r/stats/confint.html).

## Value

An object of class `beezdiscounting_linear`: a list with `subjects`
(per-unit tibble of ln k estimates, intervals and log-likelihoods), `re`
(closed-form random-effects MLEs `mu`, `sigma2`, `g`, or `NULL` if the
design is unbalanced), `design` (factor name and levels), `transform`
(boundary bookkeeping; the counts refer to the retained units), `data`
(the long frame: `id`, `x`, `y`, the factor, and the derived columns
`condition`, `unit`, `y_lin`, `d_used`, `log_jac`; `factors` may not
name one of the derived columns, and other user columns are not
carried), `conf_level`, and `call`.

## Details

Indifference points at exactly 0 or 1 are undefined under the transform;
by default they are moved to `eps` and `1 - eps` respectively
(`boundary = "clamp"`). Only exact 0 and 1 are touched — interior
points, however close to a bound, are used as observed. This is a
package decision; the paper does not address boundary values. Use
`boundary = "drop"` or `"error"` when you would rather not impute them.

## References

Hinds, D., Tegge, A. N., Stein, J. S., LaConte, S. M., McClure, S. M., &
Ferreira, M. A. R. (2026). To linearize or not to linearize: That is the
Mazur delay discounting question. *Journal of Mathematical Psychology,
130*, 103006.
[doi:10.1016/j.jmp.2026.103006](https://doi.org/10.1016/j.jmp.2026.103006)

## See also

[beezdiscounting_linear-methods](https://brentkaplan.github.io/beezdiscounting/reference/beezdiscounting_linear-methods.md)
for the S3 surface,
[`plot.beezdiscounting_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/plot.beezdiscounting_linear.md),
[`anova.beezdiscounting_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/anova.beezdiscounting_linear.md),
[`simulate_dd_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/simulate_dd_linear.md),
and
[`vignette("linearized-mazur")`](https://brentkaplan.github.io/beezdiscounting/articles/linearized-mazur.md).

## Examples

``` r
fit <- fit_dd_linear(dd_ip)
#> Warning: Clamped out-of-range y to [0, 1]: 4 value(s) > 1 set to 1, 0 value(s) < 0 set to 0.
fit
#> Linearized Mazur discounting fit (Hinds et al., 2026)
#>   100 units, 6 delays each
#>   boundary = "clamp" (eps = 0.005): 129 point(s) at 0/1, 129 clamped, 0 dropped
#>   Population ln k (mu):
#>  (all) 
#> -1.077 
#>   sigma2 = 1.090, g = 4.901, logLik(raw) = 932.52
```
