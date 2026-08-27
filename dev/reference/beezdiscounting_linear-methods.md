# Methods for linearized Mazur fits

S3 methods for objects of class `beezdiscounting_linear` returned by
[`fit_dd_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_linear.md).
All quantities are closed form; nothing is re-optimized.

## Usage

``` r
# S3 method for class 'beezdiscounting_linear'
print(x, ...)

# S3 method for class 'beezdiscounting_linear'
summary(object, ...)

# S3 method for class 'summary.beezdiscounting_linear'
print(x, ...)

# S3 method for class 'beezdiscounting_linear'
tidy(x, effects = c("subject", "population"), ...)

# S3 method for class 'beezdiscounting_linear'
glance(x, ...)

# S3 method for class 'beezdiscounting_linear'
coef(object, ...)

# S3 method for class 'beezdiscounting_linear'
confint(
  object,
  parm = c("population", "subject"),
  level = object$conf_level,
  ...
)

# S3 method for class 'beezdiscounting_linear'
augment(x, ...)

# S3 method for class 'beezdiscounting_linear'
predict(object, type = "parameters", ...)

# S3 method for class 'beezdiscounting_linear'
nobs(object, ...)

# S3 method for class 'beezdiscounting_linear'
logLik(
  object,
  level = c("population", "subject"),
  scale = c("raw", "transformed"),
  ...
)
```

## Arguments

- ...:

  Unused; present for S3 generic consistency.

- object, x:

  A `beezdiscounting_linear` fit (for `print.summary()`, the object
  returned by [`summary()`](https://rdrr.io/r/base/summary.html)).

- effects:

  For [`tidy()`](https://generics.r-lib.org/reference/tidy.html):
  `"subject"` (default) returns the per-unit table; `"population"`
  returns the random-effects MLEs `mu_<level>`, `sigma2` and `g`.

- parm:

  Which interval [`confint()`](https://rdrr.io/r/stats/confint.html)
  returns: `"population"` (default) is an ANOVA-style t-interval on each
  condition mean, using the between-unit mean square with `N - C`
  degrees of freedom; `"subject"` gives the per-subject t-intervals on
  ln k. It does not select parameter names.

- level:

  For [`confint()`](https://rdrr.io/r/stats/confint.html), the
  confidence level; defaults to the fit's `conf_level`, so
  `confint(parm = "subject")` reproduces the intervals in
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  `fit$subjects` unless a different `level` is given. For
  [`logLik()`](https://rdrr.io/r/stats/logLik.html), which likelihood to
  return: `"population"` (the random-effects MLE, `df` = C + 2) or
  `"subject"` (the sum of the per-unit log-likelihoods, `df` = 2 per
  unit).

- type:

  For [`predict()`](https://rdrr.io/r/stats/predict.html): only
  `"parameters"` exists, returning the per-subject `ln k` and `k` with
  their t-intervals (`id`, `condition`, `logk`, `logk_lower`,
  `logk_upper`, `k`, `k_lower`, `k_upper`). There is no `newdata`
  prediction;
  [`augment()`](https://generics.r-lib.org/reference/augment.html) gives
  the fitted curve.

- scale:

  For [`logLik()`](https://rdrr.io/r/stats/logLik.html): `"raw"`
  (default) is Jacobian-corrected back to the indifference-point scale;
  `"transformed"` is on the linearized scale.

## Value

[`print()`](https://rdrr.io/r/base/print.html) and `print.summary()`
return their input invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.beezdiscounting_linear`;
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html),
[`augment()`](https://generics.r-lib.org/reference/augment.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) return tibbles.
[`augment()`](https://generics.r-lib.org/reference/augment.html) adds
`.fitted` and `.resid` on the raw indifference-point scale and
`.std_resid`, the transformed-scale residual `y_lin - ln k_i` divided by
the model's error standard deviation (`sqrt(sigma2)` from the
random-effects fit, or the pooled within-subject mean square when that
component is unavailable); `.std_resid` is `NA` for points dropped under
`boundary = "drop"`. [`coef()`](https://rdrr.io/r/stats/coef.html)
returns the named vector of condition means `mu`;
[`confint()`](https://rdrr.io/r/stats/confint.html) returns a two-column
matrix of lower/upper bounds;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
usable transformed observations;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns a `"logLik"`
object with `df` and `nobs` attributes. Methods needing the
random-effects component error when the design was unbalanced and `re`
is `NULL`.

## Details

[`confint()`](https://rdrr.io/r/stats/confint.html) uses `parm` to
SELECT WHICH INTERVAL is returned, not to filter parameter names as in
[`stats::confint()`](https://rdrr.io/r/stats/confint.html).
[`logLik()`](https://rdrr.io/r/stats/logLik.html) reports the raw-scale
(Jacobian-corrected) log-likelihood by default, which makes it
comparable with this package's Gaussian NLS and TMB log-likelihoods.

## See also

[`fit_dd_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_linear.md),
[`anova.beezdiscounting_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/anova.beezdiscounting_linear.md),
[`plot.beezdiscounting_linear()`](https://brentkaplan.github.io/beezdiscounting/reference/plot.beezdiscounting_linear.md)
