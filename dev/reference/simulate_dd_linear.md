# Simulate indifference points from the linearized Mazur random-effects model

Draws `theta_ic ~ N(mu_c, g sigma2 / T)`, `y_ijc ~ N(theta_ic, sigma2)`
and returns `D = 1 / (1 + exp(y) t)` (Hinds et al., 2026, Sec. 4.1),
which is always strictly inside (0, 1).

## Usage

``` r
simulate_dd_linear(
  n_subjects,
  delays,
  mu,
  sigma2,
  g,
  seed = NULL,
  attach_truth = FALSE
)
```

## Arguments

- n_subjects:

  Integer, scalar (recycled) or one per condition.

- delays:

  Positive numeric vector of delays (`T = length(delays)`).

- mu:

  Numeric vector of population mean ln k per condition; names become
  condition labels.

- sigma2:

  Measurement-error variance on the transformed scale.

- g:

  Random-effect variance multiplier: `Var(theta) = g * sigma2 / T`.

- seed:

  Optional seed; the global RNG state is restored afterwards.

- attach_truth:

  If `TRUE`, attach `attr(, "truth")`.

## Value

Long `data.frame(id, condition, x, y)`.

## Examples

``` r
s <- simulate_dd_linear(30, c(7, 30, 180, 365), mu = c(A = -6, B = -5),
                        sigma2 = 2, g = 10, seed = 1)
anova(fit_dd_linear(s, factors = "condition"))
#> # A tibble: 1 × 6
#>   hypothesis     F   df1   df2 p_value cohens_d
#>   <chr>      <dbl> <int> <dbl>   <dbl>    <dbl>
#> 1 A = B       4.59     1    58  0.0363   -0.563
```
