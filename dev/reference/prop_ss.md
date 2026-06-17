# Calculate proportion of SIR/SS responses at each k value

Calculate proportion of SIR/SS responses at each k value

## Usage

``` r
prop_ss(dat)
```

## Arguments

- dat:

  Dataframe (longform) with subjectid, questionid, and response (0 for
  SIR/SS and 1 for LDR/LL)

## Value

Dataframe with proportion of SIR/SS responses at each k rank

## Examples

``` r
prop_ss(mcq27)
#> # A tibble: 9 × 2
#>    k_rank prop_ss
#>     <dbl>   <dbl>
#> 1 0.00016    1   
#> 2 0.0004     0.83
#> 3 0.001      0.67
#> 4 0.0025     0.5 
#> 5 0.006      0.5 
#> 6 0.016      0.5 
#> 7 0.041      0.33
#> 8 0.1        0   
#> 9 0.25       0   
```
