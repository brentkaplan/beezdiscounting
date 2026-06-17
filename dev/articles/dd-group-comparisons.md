# Comparing discounting rates between groups

A common discounting question is not “what is k?” but “does k differ
between groups?” – treatment vs. control, clinical vs. community,
commodity A vs. commodity B. Both modeling tiers in beezdiscounting
answer it the same way: the discount rate is estimated on the
natural-log scale, so a between-subject design enters as a linear model
on log k (`log k = X beta`), and every marginal mean or contrast is a
linear combination of coefficients.
[`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)
and
[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)
work on top of that shared backbone for:

- [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)
  – mixed-effects indifference-point models (Wald inference),
- [`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md)
  – the structural trial-level choice model (Wald inference),
- [`fit_dd_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_brms.md)
  and
  [`fit_dd_choice_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice_brms.md)
  – their Bayesian counterparts (posterior draws; see the last section).

The TMB examples below run at build time; the Bayesian section shows
code without running it (Stan sampling takes minutes – copy the chunks
into your session).

## Simulate a two-group study

[`simulate_dd_ip()`](https://brentkaplan.github.io/beezdiscounting/reference/simulate_dd_ip.md)
generates indifference points with a known between-condition shift in
log k. Here condition `C2`’s k is three times `C1`’s
(`delta_k = log(3)`). To have a second factor for the grouped contrasts
shown later, we also assign each subject a (null) `site`.

``` r

library(beezdiscounting)

sim <- simulate_dd_ip(
  n_subjects = 24, n_conditions = 2, delta_k = c(0, log(3)), seed = 42
)
# conditions alternate across subjects, so assign site in blocks of two to
# get the full condition x site crossing
subjects <- unique(sim$id)
sim$site <- factor(ifelse(
  ((match(sim$id, subjects) - 1) %/% 2) %% 2 == 0, "siteA", "siteB"
))
with(unique(sim[, c("id", "condition", "site")]), table(condition, site))
#>          site
#> condition siteA siteB
#>        C1     6     6
#>        C2     6     6
head(sim)
#> # A tibble: 6 × 5
#>   id    condition     x      y site 
#>   <fct> <fct>     <dbl>  <dbl> <fct>
#> 1 1     C1            7 0.990  siteA
#> 2 1     C1           30 0.648  siteA
#> 3 1     C1          180 0.127  siteA
#> 4 1     C1          365 0.0523 siteA
#> 5 1     C1          730 0.0194 siteA
#> 6 1     C1         1460 0.0429 siteA
```

## Fit with a factor design on log k

`factors` puts the design on log k; every other parameter stays
population-level. (`multi_start = FALSE` keeps the vignette quick; leave
the default on for real analyses.)

``` r

fit <- fit_dd_tmb(
  sim,
  factors = c("condition", "site"),
  multi_start = FALSE, verbose = 0
)
tidy(fit)
#> # A tibble: 5 × 8
#>   term           estimate std.error statistic   p.value component estimate_scale
#>   <chr>             <dbl>     <dbl>     <dbl>     <dbl> <chr>     <chr>         
#> 1 k:(Intercept)   0.00923   0.00253    -17.1   1.09e-65 fixed     natural       
#> 2 k:conditionC2   2.97      0.813        3.98  6.90e- 5 fixed     natural       
#> 3 k:sitesiteB     1.47      0.403        1.41  1.57e- 1 fixed     natural       
#> 4 sigma_u (log1…  0.312    NA           NA    NA        variance  log10         
#> 5 phi (precisio… 10.5      NA           NA    NA        variance  natural       
#> # ℹ 1 more variable: term_display <chr>
```

## Estimated marginal means of k

[`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)
averages the fitted design over the factor levels you are not looking at
(equal weights, the emmeans convention) and back-transforms:
`k = exp(k_log)`, with the standard error and the Wald interval formed
on the log scale.

``` r

get_dd_param_emms(fit)
#> # A tibble: 4 × 6
#>   level                          k k_log std.error conf.low conf.high
#>   <chr>                      <dbl> <dbl>     <dbl>    <dbl>     <dbl>
#> 1 condition=C1, site=siteA 0.00923 -4.68     0.274  0.00540    0.0158
#> 2 condition=C1, site=siteB 0.0136  -4.30     0.276  0.00791    0.0234
#> 3 condition=C2, site=siteA 0.0274  -3.60     0.277  0.0159     0.0472
#> 4 condition=C2, site=siteB 0.0404  -3.21     0.280  0.0233     0.0700
```

Marginalize over `site` to get per-condition means, or condition on a
level with `at`:

``` r

get_dd_param_emms(fit, factors_in_emm = "condition")
#> # A tibble: 2 × 6
#>   level             k k_log std.error conf.low conf.high
#>   <chr>         <dbl> <dbl>     <dbl>    <dbl>     <dbl>
#> 1 condition=C1 0.0112 -4.49     0.224  0.00722    0.0174
#> 2 condition=C2 0.0333 -3.40     0.229  0.0213     0.0521

get_dd_param_emms(fit, at = list(site = "siteA"))
#> # A tibble: 2 × 6
#>   level                          k k_log std.error conf.low conf.high
#>   <chr>                      <dbl> <dbl>     <dbl>    <dbl>     <dbl>
#> 1 condition=C1, site=siteA 0.00923 -4.68     0.274  0.00540    0.0158
#> 2 condition=C2, site=siteA 0.0274  -3.60     0.277  0.0159     0.0472
```

## Contrasts

[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)
turns the same reference grid into factor-level contrasts. Estimates are
reported on the log10 scale (the customary discounting metric) and,
because differences in log k are log ratios, also as multiplicative
ratios: a ratio of 3 means one group discounts three times as steeply.
P-values are adjusted with
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) (`"holm"`
by default).

``` r

cmp <- get_dd_comparisons(fit, compare_specs = ~condition)
cmp$k$contrasts_log10
#> # A tibble: 1 × 8
#>   contrast         estimate std.error statistic    df conf.low conf.high p.value
#>   <chr>               <dbl>     <dbl>     <dbl> <dbl>    <dbl>     <dbl>   <dbl>
#> 1 condition=C1 - …   -0.473     0.139     -3.41   Inf   -0.745    -0.201 6.46e-4
cmp$k$contrasts_ratio
#> # A tibble: 1 × 5
#>   contrast                    ratio conf.low conf.high  p.value
#>   <chr>                       <dbl>    <dbl>     <dbl>    <dbl>
#> 1 condition=C1 - condition=C2 0.336    0.180     0.629 0.000646
```

The simulated truth was `k(C2) = 3 k(C1)`, so the displayed `C1 - C2`
ratio should land near 1/3 – and it does. With more than one factor,
`contrast_by` conditions the contrasts on another factor’s levels
(adjustment is applied per by-cell; with no interaction in the model the
per-site contrasts coincide, as here):

``` r

cmp_by <- get_dd_comparisons(
  fit,
  compare_specs = ~ condition + site, contrast_by = "site"
)
cmp_by$k$contrasts_log10
#> # A tibble: 2 × 9
#>   site  contrast   estimate std.error statistic    df conf.low conf.high p.value
#>   <chr> <chr>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>     <dbl>   <dbl>
#> 1 siteA condition…   -0.473     0.139     -3.41   Inf   -0.745    -0.201 6.46e-4
#> 2 siteB condition…   -0.473     0.139     -3.41   Inf   -0.745    -0.201 6.46e-4
```

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) flattens any
comparison object into one cross-backend frame:

``` r

tidy(cmp)
#> # A tibble: 1 × 9
#>   param contrast   estimate std.error statistic    df conf.low conf.high p.value
#>   <chr> <chr>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>     <dbl>   <dbl>
#> 1 k     condition…   -0.473     0.139     -3.41   Inf   -0.745    -0.201 6.46e-4
```

## Trial-level choice models

The structural choice model estimates the same log-k design from binary
smaller-sooner vs. larger-later choices, so the comparison surface is
identical. (The simulator has no built-in group effect; stack two
simulated groups.)

``` r

ctrl <- simulate_dd_choice(n_subjects = 12, log_k_pop = log(0.01), seed = 7)
trt <- simulate_dd_choice(n_subjects = 12, log_k_pop = log(0.04), seed = 8)
trt$id <- paste0("t", trt$id)
ctrl$group <- "ctrl"
trt$group <- "treat"
choice_data <- rbind(ctrl, trt)

fit_choice <- fit_dd_choice(
  choice_data,
  mode = "structural", factors = "group", verbose = 0
)
get_dd_param_emms(fit_choice)
#> # A tibble: 2 × 6
#>   level            k k_log std.error conf.low conf.high
#>   <chr>        <dbl> <dbl>     <dbl>    <dbl>     <dbl>
#> 1 group=ctrl  0.0113 -4.48     0.234  0.00716    0.0179
#> 2 group=treat 0.0369 -3.30     0.212  0.0243     0.0559
get_dd_comparisons(fit_choice)$k$contrasts_ratio
#> # A tibble: 1 × 5
#>   contrast                 ratio conf.low conf.high  p.value
#>   <chr>                    <dbl>    <dbl>     <dbl>    <dbl>
#> 1 group=ctrl - group=treat 0.307    0.167     0.566 0.000152
```

## The Bayesian analogs

[`fit_dd_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_brms.md)
and
[`fit_dd_choice_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice_brms.md)
take the same `factors`/`factor_interaction`/`continuous_covariates`
arguments, and
[`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)/[`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)
accept the fits directly:

``` r

fit_b <- fit_dd_brms(
  sim,
  equation = "mazur", factors = c("condition", "site"),
  chains = 4, cores = 4, seed = 1
)
get_dd_param_emms(fit_b, factors_in_emm = "condition")
get_dd_comparisons(fit_b, compare_specs = ~condition)

fit_cb <- fit_dd_choice_brms(
  choice_data,
  equation = "mazur", factors = "group",
  chains = 4, cores = 4, seed = 1
)
get_dd_comparisons(fit_cb)
```

The tables have the same shape, with draws-based summaries in place of
Wald machinery:

- Estimates are posterior medians; intervals are equal-tailed quantiles
  of the (transformed) draws – the marginal means and contrasts are
  computed per posterior draw over the same reference grid, then
  summarized, so back-transformation is exact (no delta method).
- `statistic`, `df`, and `p.value` are `NA`. In their place `post.prob`
  reports the posterior probability of direction: the probability that
  the contrast is positive (or negative, whichever is larger). A
  `post.prob` near 1 plays the practical role of a small p-value.
- There is no multiplicity adjustment, by design: the contrasts are
  deterministic functions of one joint posterior, which already encodes
  their dependence. Asking for `adjust` produces a warning rather than
  an adjusted column.

## Which backend when

Use the TMB fitters for fast iteration and conventional Wald tests with
familiar adjusted p-values. Reach for the brms fitters when you want
full posterior uncertainty for derived quantities (the EMM and ratio
intervals above need no asymptotics), prior information about plausible
discounting rates, or small-sample inference you can defend without
appealing to asymptotic normality. The vignette
[`vignette("bayesian-discounting")`](https://brentkaplan.github.io/beezdiscounting/articles/bayesian-discounting.md)
covers the Bayesian tier itself – priors, families, convergence – in
more depth.
