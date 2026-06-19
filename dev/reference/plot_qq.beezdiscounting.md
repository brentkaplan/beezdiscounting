# Random-effect normal QQ plots for discounting models

QQ plot of the estimated subject random-effect deviates against a normal
reference – the standard check on the Gaussian random-effects
assumption. Methods exist for the TMB indifference-point
([`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md))
and choice
([`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md))
models. Bayesian fits are intentionally excluded; use
[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
and MCMC diagnostics there instead. The plotted deviates are shrunken
(empirical-Bayes) estimates, so the normal reference is approximate.

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
plot_qq(object, which = NULL, ...)

# S3 method for class 'beezdiscounting_choice'
plot_qq(object, which = NULL, ...)
```

## Arguments

- object:

  A fitted `beezdiscounting_tmb` or `beezdiscounting_choice` model.

- which:

  Optional character vector selecting random-effect terms to show.

- ...:

  Unused.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
