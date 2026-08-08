# Calculate proportion of SIR/SS responses at each k value

Calculate proportion of SIR/SS responses at each k value

## Usage

``` r
prop_ss(dat, items = 27)
```

## Arguments

- dat:

  Dataframe (longform) with subjectid, questionid, and response (0 for
  SIR/SS and 1 for LDR/LL)

- items:

  Number of MCQ items (27 or 21)

## Value

Dataframe with proportion of SIR/SS responses at each k rank

## Details

`items` must match the instrument actually administered. Question ids
1-21 are valid in both the 21- and 27-item designs, so passing the wrong
`items` does not error – it silently pools responses into the wrong
k-rank rows. If the observed question ids do not exactly match the
requested design, `prop_ss()` warns.

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
dat21 <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
prop_ss(dat21, items = 21)
#> # A tibble: 7 × 2
#>   k_rank prop_ss
#>    <dbl>   <dbl>
#> 1      1       0
#> 2      2       0
#> 3      3       0
#> 4      4       0
#> 5      5       0
#> 6      6       0
#> 7      7       0
```
