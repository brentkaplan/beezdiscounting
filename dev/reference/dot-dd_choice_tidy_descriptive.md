# Tidy a descriptive (Young 2018) choice model into a coefficient tibble

Fixed-effect (`theta`) rows are the logit-scale sensitivities and are
ALWAYS on the identity scale - never exponentiated across report spaces
(mirroring the structural `beta0`). Variance rows report the
random-slope SDs and their correlation on the natural scale. Wald
`statistic`/`p.value` are computed on the estimation scale.

## Usage

``` r
.dd_choice_tidy_descriptive(x, effects, ...)
```

## Arguments

- x:

  A `beezdiscounting_choice` object (descriptive mode).

- effects:

  Resolved character vector (`"fixed"` / `"ran_pars"`).

- ...:

  Unused.

## Value

A tibble with the standard 8-column broom contract.
