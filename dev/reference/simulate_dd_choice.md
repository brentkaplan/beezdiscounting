# Simulate trial-level SS-vs-LL choice data

Generates binary choices from either the structural or descriptive
choice model. When `mode = "structural"` (default) choices are drawn
from the structural model:
`P(LL) = plogis(beta0 + gamma * ((ll/ss) * D(k, delay) - 1))`,
`k_i = exp(log_k_pop + sigma_u * u_i)`. When `mode = "descriptive"`
choices follow Young's (2018) correlated random-slope logistic model:
each subject receives per-subject slopes drawn from a bivariate normal
with mean `theta` and covariance `Sigma`, and
`P(LL) = plogis((theta[1] + b_i[1]) * log(ll/ss) + (theta[2] + b_i[2]) * log(delay + 1))`.
The descriptive branch is the known-truth generator used by recovery
tests.

## Usage

``` r
simulate_dd_choice(
  n_subjects = 50,
  ss_amount = c(40, 55, 31, 14, 47, 25, 78, 40, 11, 67),
  ll_amount = c(65, 75, 85, 25, 60, 30, 80, 65, 30, 75),
  delay = c(25, 61, 14, 19, 160, 7, 4, 62, 7, 119),
  mode = c("structural", "descriptive"),
  log_k_pop = log(0.02),
  sigma_u = 0.5,
  gamma = 4,
  beta0 = 0,
  equation = c("mazur", "exponential"),
  theta = c(1.5, -0.4),
  re_sd = c(0.5, 0.3),
  re_cor = -0.2,
  return_truth = FALSE,
  seed = NULL
)
```

## Arguments

- n_subjects:

  Integer number of subjects.

- ss_amount, ll_amount, delay:

  Numeric trial-design vectors (same length); each subject is presented
  all trials. Defaults are a built-in grid.

- mode:

  `"structural"` (default) or `"descriptive"`. Selects the generative
  model.

- log_k_pop, sigma_u, gamma, beta0:

  Structural truth (used only when `mode = "structural"`).

- equation:

  `"mazur"` or `"exponential"` (structural mode only).

- theta:

  Length-2 numeric vector of population-level slopes for the descriptive
  model: `theta[1]` on `log(ll/ss)`, `theta[2]` on `log(delay + 1)`.

- re_sd:

  Length-2 numeric vector of random-effect standard deviations
  (descriptive mode only).

- re_cor:

  Correlation between the two random slopes (descriptive mode only).

- return_truth:

  Logical. If `TRUE` and `mode = "descriptive"`, the returned tibble
  carries two attributes: `subject_slopes` (n_subjects x 2 matrix of
  realized `b_i` values) and `Sigma` (the 2x2 covariance matrix).

- seed:

  Optional integer seed for reproducibility.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with
columns `id`, `ss_amount`, `ll_amount`, `delay`, `choice` (0/1, 1 = LL
chosen). When `mode = "descriptive"` and `return_truth = TRUE`, the
tibble also carries `subject_slopes` and `Sigma` attributes.

## Examples

``` r
# Structural (default)
sim <- simulate_dd_choice(n_subjects = 20, seed = 1)
head(sim)
#> # A tibble: 6 × 5
#>   id    ss_amount ll_amount delay choice
#>   <chr>     <dbl>     <dbl> <dbl>  <int>
#> 1 1            40        65    25      0
#> 2 1            55        75    61      0
#> 3 1            31        85    14      1
#> 4 1            14        25    19      1
#> 5 1            47        60   160      0
#> 6 1            25        30     7      0

# Descriptive (Young 2018 correlated random-slope model)
sim_desc <- simulate_dd_choice(n_subjects = 20, mode = "descriptive",
                               theta = c(1.5, -0.4), seed = 42)
head(sim_desc)
#> # A tibble: 6 × 5
#>   id    ss_amount ll_amount delay choice
#>   <chr>     <dbl>     <dbl> <dbl>  <int>
#> 1 1            40        65    25      0
#> 2 1            55        75    61      0
#> 3 1            31        85    14      0
#> 4 1            14        25    19      0
#> 5 1            47        60   160      0
#> 6 1            25        30     7      0
```
