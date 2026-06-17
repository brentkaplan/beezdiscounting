# Summarize a TMB mixed-effects discounting fit

Summarize a TMB mixed-effects discounting fit

## Usage

``` r
# S3 method for class 'beezdiscounting_tmb'
summary(object, report_space = c("natural", "log10", "internal", "log"), ...)
```

## Arguments

- object:

  A `beezdiscounting_tmb` object.

- report_space:

  Scale for the fixed-effect (`beta_k`) estimates and standard errors in
  the coefficient table: `"natural"` (default; `k` on the natural scale
  via [`exp()`](https://rdrr.io/r/base/Log.html)), `"log10"` (log10-k),
  `"internal"` or `"log"` (log-k, the estimation scale; the two coincide
  for `beta_k`). `statistic` and `p.value` are always computed on the
  estimation (log-k) scale regardless of `report_space`.

- ...:

  Unused.

## Value

An object of class `"summary.beezdiscounting_tmb"` with components:
`call`, `model_class`, `backend`, `equation`, `family`, `coefficients`,
`variance_components`, `n_subjects`, `nobs`, `converged`, `logLik`,
`AIC`, `BIC`, `notes`.

## Examples

``` r
# \donttest{
dd <- simulate_dd_ip(n_subjects = 20, seed = 1)
fit <- fit_dd_tmb(dd, equation = "mazur", family = "sltb",
                  random_effects = k ~ 1, verbose = 0)
summary(fit)
#> 
#> TMB Mixed-Effects Discounting Model Summary
#> ================================================== 
#> 
#> Call:
#> fit_dd_tmb(data = dd, equation = "mazur", family = "sltb", random_effects = k ~ 
#>     1, verbose = 0) 
#> 
#> Equation: mazur 
#> Family: sltb 
#> Backend: TMB_mixed 
#> Convergence: Yes 
#> Subjects: 20  Observations: 140 
#> 
#> --- Fixed Effects (k) ---
#>           term estimate std.error statistic p.value
#>  k:(Intercept)   0.0108    0.0015  -33.6239  <2e-16
#> 
#> --- Variance Components ---
#>                Component Estimate   Scale
#>  sigma_u (log10-k RE SD)   0.2309   log10
#>          phi (precision)  12.4798 natural
#> 
#> --- Fit Statistics ---
#> Log-likelihood: 198.92 
#> AIC: -391.84   BIC: -383.01 
summary(fit, report_space = "log10")
#> 
#> TMB Mixed-Effects Discounting Model Summary
#> ================================================== 
#> 
#> Call:
#> fit_dd_tmb(data = dd, equation = "mazur", family = "sltb", random_effects = k ~ 
#>     1, verbose = 0) 
#> 
#> Equation: mazur 
#> Family: sltb 
#> Backend: TMB_mixed 
#> Convergence: Yes 
#> Subjects: 20  Observations: 140 
#> 
#> --- Fixed Effects (log10 k) ---
#>           term estimate std.error statistic p.value
#>  k:(Intercept)  -1.9682    0.0585  -33.6239  <2e-16
#> 
#> --- Variance Components ---
#>                Component Estimate   Scale
#>  sigma_u (log10-k RE SD)   0.2309   log10
#>          phi (precision)  12.4798 natural
#> 
#> --- Fit Statistics ---
#> Log-likelihood: 198.92 
#> AIC: -391.84   BIC: -383.01 
# }
```
