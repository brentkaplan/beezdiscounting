# Calculate scores, answers, and timing for 5.5 trial probability discounting from Qualtrics template

Calculate scores, answers, and timing for 5.5 trial probability
discounting from Qualtrics template

## Usage

``` r
calc_pd(df)
```

## Arguments

- df:

  A dataframe containing all the columns from the template.

## Value

A dataframe with h/ep50 values, answers, timing

## Examples

``` r
calc_pd(five.fivetrial_pd)
#> # A tibble: 21 × 12
#>    ResponseId index     q firstclick lastclick pagesubmit totalclicks response
#>         <int> <chr> <dbl>      <dbl>     <dbl>      <dbl>       <dbl> <chr>   
#>  1          1 I16       1       3.98      3.98       5.18           1 sc      
#>  2          1 I8        2       4.01      4.01       4.76           1 lu      
#>  3          1 I12       3       2.06      2.06       3.25           1 sc      
#>  4          1 I10       4       1.52      1.52       3.02           1 sc      
#>  5          1 I9        5       2.25      2.95       3.74           2 lu      
#>  6          2 I16       1       2.87      2.87       3.88           1 sc      
#>  7          2 I8        2       3.74      3.74       4.86           1 sc      
#>  8          2 I4        3       1.16      1.16       6.36           1 sc      
#>  9          2 I2        4       3.06      3.06       5.41           1 sc      
#> 10          2 I1        5       2.05      2.05       5.10           1 sc      
#> # ℹ 11 more rows
#> # ℹ 4 more variables: attentionflag <chr>, hval <dbl>, etheta50 <dbl>,
#> #   ep50 <dbl>
```
