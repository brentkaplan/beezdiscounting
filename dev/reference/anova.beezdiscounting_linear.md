# F-test for condition means in a linearized Mazur fit (Hinds et al. 2026, Prop. 3.5)

The paper's F-test is exact when the random-effects variance estimate
`g-hat > 0` (Prop. 3.5's assumption); when `g-hat = 0` the statistic is
still reported, with a warning.

## Usage

``` r
# S3 method for class 'beezdiscounting_linear'
anova(object, hypothesis = NULL, pairwise = FALSE, ...)
```

## Arguments

- object:

  A `beezdiscounting_linear` fit with a factor.

- hypothesis:

  `NULL` (all condition means equal) or a list of character vectors,
  each naming levels constrained equal under H0.

- pairwise:

  If `TRUE`, test every pair of levels (uncorrected p-values); cannot be
  combined with `hypothesis`.

- ...:

  Unused; present for S3 generic consistency.

## Value

Tibble of class `beezdiscounting_linear_anova`.

## See also

[beezdiscounting_linear-methods](https://brentkaplan.github.io/beezdiscounting/reference/beezdiscounting_linear-methods.md)
for the other S3 methods.
