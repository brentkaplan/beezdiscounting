# Factor-level comparisons of the discount rate `k`

Computes factor-level contrasts of the discount rate `k` from a fitted
`beezdiscounting_tmb` model (or a structural `beezdiscounting_choice`
fit). `k` is linear in the fixed-effect coefficients on the natural-log
scale (`log k = X beta_k`), so each contrast is a linear combination of
`beta_k` with a Wald standard error from the `beta_k` block of
[`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html)'s
fixed-effect covariance. Contrasts are **reported on the `log10` scale**
and, optionally, as multiplicative ratios (`ratio = exp(est_log)`). The
returned container mirrors the beezdemand `beezdemand_comparison` shape,
so
[`tidy.beezdiscounting_comparison()`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_comparison.md)
gives a flat, cross-backend frame. Bayesian fits
(`beezdiscounting_brms`, `beezdiscounting_choice_brms`) compute per-draw
contrasts over the same reference grid: posterior medians with quantile
credible intervals and `post.prob` in place of adjusted p-values (no
multiplicity adjustment; the joint posterior already encodes contrast
dependence).

## Usage

``` r
get_dd_comparisons(
  fit,
  compare_specs = NULL,
  contrast_type = c("pairwise", "trt.vs.ctrl"),
  contrast_by = NULL,
  adjust = "holm",
  at = NULL,
  ci_level = 0.95,
  report_ratios = TRUE,
  ...
)
```

## Arguments

- fit:

  A `beezdiscounting_tmb`, structural `beezdiscounting_choice`,
  `beezdiscounting_brms`, or `beezdiscounting_choice_brms` object.

- compare_specs:

  Optional one-sided formula naming the factor subset to contrast (e.g.
  `~ condition`). Omitted fitted factors are marginalized over with
  equal weights across the full crossing of their levels. If `NULL`
  (default), all fitted factors are retained. Unknown names abort.

- contrast_type:

  Character. `"pairwise"` (all `choose(n, 2)` pairs, factor-level order)
  or `"trt.vs.ctrl"` (each level vs. the first/reference level).

- contrast_by:

  Optional `NULL` (default) or character vector of factor name(s) within
  `compare_specs` to condition the contrasts on. Within each observed
  combination of by-level(s), pairwise (or `trt.vs.ctrl`) contrasts are
  computed over the remaining (non-by) factors, with p-value adjustment
  applied **per by-cell**. A `contrast_by` factor absent from
  `compare_specs` aborts.

- adjust:

  Character. P-value adjustment method; must be one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html)
  (default `"holm"`). emmeans-only methods such as `"tukey"`/`"sidak"`
  are rejected (the TMB backend uses an asymptotic *z* +
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)).

- at:

  Named list specifying factor levels and/or continuous-covariate values
  to condition on, as in
  [`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md).

- ci_level:

  Numeric. Confidence level for intervals (default 0.95).

- report_ratios:

  Logical. If `TRUE` (default), include a `contrasts_ratio` block
  (multiplicative ratios).

- ...:

  Additional arguments (reserved; `factors_in_emm` is accepted as a
  lower-level alternative to `compare_specs`).

## Value

A `beezdiscounting_comparison` object: a list with a single element `k`,
itself a list with `emmeans` (the EMM tibble from
[`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)),
`contrasts_log10` (columns `contrast`, `estimate`, `std.error`,
`statistic`, `df`, `conf.low`, `conf.high`, `p.value`), and (if
`report_ratios`) `contrasts_ratio` (columns `contrast`, `ratio`,
`conf.low`, `conf.high`, `p.value`). When `contrast_by` is active, the
contrast tables gain leading by-column(s) before `contrast`. Attributes
`backend`, `adjustment_method`, `compare_specs_used`,
`contrast_type_used`, and `contrast_by_used` describe the call.

## See also

[`tidy.beezdiscounting_comparison()`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_comparison.md)
for the cross-backend frame.

## Examples

``` r
# \donttest{
sim <- simulate_dd_ip(
  n_subjects = 30, n_conditions = 2, delta_k = c(0, log(3)), seed = 1
)
fit <- fit_dd_tmb(sim, factors = "condition", verbose = 0)
get_dd_comparisons(fit, contrast_type = "pairwise")
#> $k
#> $k$emmeans
#> # A tibble: 2 × 6
#>   level             k k_log std.error conf.low conf.high
#>   <chr>         <dbl> <dbl>     <dbl>    <dbl>     <dbl>
#> 1 condition=C1 0.0112 -4.49     0.152  0.00834    0.0151
#> 2 condition=C2 0.0276 -3.59     0.154  0.0204     0.0373
#> 
#> $k$contrasts_log10
#> # A tibble: 1 × 8
#>   contrast         estimate std.error statistic    df conf.low conf.high p.value
#>   <chr>               <dbl>     <dbl>     <dbl> <dbl>    <dbl>     <dbl>   <dbl>
#> 1 condition=C1 - …   -0.391    0.0935     -4.18   Inf   -0.574    -0.207 2.95e-5
#> 
#> $k$contrasts_ratio
#> # A tibble: 1 × 5
#>   contrast                    ratio conf.low conf.high   p.value
#>   <chr>                       <dbl>    <dbl>     <dbl>     <dbl>
#> 1 condition=C1 - condition=C2 0.407    0.267     0.620 0.0000295
#> 
#> 
#> attr(,"class")
#> [1] "beezdiscounting_comparison"
#> attr(,"backend")
#> [1] "tmb"
#> attr(,"adjustment_method")
#> [1] "holm"
#> attr(,"compare_specs_used")
#> [1] "all fitted factors"
#> attr(,"contrast_type_used")
#> [1] "pairwise"
#> attr(,"contrast_by_used")
#> [1] "NULL"
# }
```
