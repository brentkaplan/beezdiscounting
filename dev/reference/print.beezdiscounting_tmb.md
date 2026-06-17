# Print a TMB mixed-effects discounting fit

Print a TMB mixed-effects discounting fit

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
print(x, ...)
```

## Arguments

- x:

  A `beezdiscounting_tmb` object.

- ...:

  Unused.

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
dd <- simulate_dd_ip(n_subjects = 20, seed = 1)
fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb",
                  random_effects = k ~ 1, verbose = 0)
print(fit)
#> 
#> TMB Mixed-Effects Discounting Model
#> 
#> Call:
#> fit_dd_tmb(data = dd, equation = "mazur", family = "sltb", random_effects = k ~ 
#>     1, verbose = 0)
#> 
#> Equation: mazur 
#> Family: sltb 
#> Convergence: Yes 
#> Number of subjects: 20 
#> Number of observations: 140 
#> Random effects: 1 (k ~ 1)
#> Log-likelihood: 198.92 
#> AIC: -391.84 
#> 
#> Fixed Effects (log k):
#> k:(Intercept) 
#>        -4.532 
#> 
#> Use summary() for full results.
# }
```
