# Provide a summary of the results from the MCQ output table.

Provide a summary of the results from the MCQ output table.

## Usage

``` r
summarize_mcq(res, na.rm = TRUE)
```

## Arguments

- res:

  Dataframe with MCQ results (output from the `calc_mcq` function)

- na.rm:

  Boolean whether to remove NAs from the calculation

## Value

Dataframe with summary statistics

## Examples

``` r
summarize_mcq(score_mcq27(mcq27))
#> # A tibble: 10 × 4
#>    Metric                  Mean     SD    SEM
#>    <chr>                  <dbl>  <dbl>  <dbl>
#>  1 overall_k             0.0328 0.0458 0.0324
#>  2 small_k               0.0131 0.0176 0.0125
#>  3 medium_k              0.0326 0.0439 0.0311
#>  4 large_k               0.0326 0.0457 0.0323
#>  5 geomean_k             0.0240 0.0330 0.0233
#>  6 overall_consistency   0.963  0      0     
#>  7 small_consistency     1      0      0     
#>  8 medium_consistency    1      0      0     
#>  9 large_consistency     1      0      0     
#> 10 composite_consistency 1      0      0     
```
