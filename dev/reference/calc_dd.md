# Calculate scores, answers, and timing for 5.5 trial delay discounting from Qualtrics template

Calculate scores, answers, and timing for 5.5 trial delay discounting
from Qualtrics template

## Usage

``` r
calc_dd(df)
```

## Arguments

- df:

  A dataframe containing all the columns from the template.

## Value

A dataframe with k/ed50 values, answers, timing

## Examples

``` r
calc_dd(five.fivetrial_dd)
#> # A tibble: 22 × 11
#>    ResponseId index     q firstclick lastclick pagesubmit totalclicks response
#>         <int> <chr> <dbl>      <dbl>     <dbl>      <dbl>       <dbl> <chr>   
#>  1          1 I16       1      1.76      1.76        3.34           1 ll      
#>  2          1 I24       2      7.73      7.73        8.46           1 ss      
#>  3          1 I20       3      1.56      1.56        3.38           1 ll      
#>  4          1 I22       4      2.33      3.95        4.50           2 ss      
#>  5          1 I21       5      3.16      3.16        3.73           1 ss      
#>  6          2 I16       1      3.78      3.78        4.35           1 ss      
#>  7          2 I8        2      1.45      1.45        3.19           1 ss      
#>  8          2 I4        3      1.18      1.18        3.14           1 ll      
#>  9          2 I6        4      0.873     0.873       3.26           1 ss      
#> 10          2 I5        5      2.62      2.62        3.26           1 ss      
#> # ℹ 12 more rows
#> # ℹ 3 more variables: attentionflag <chr>, kval <dbl>, ed50 <dbl>
```
