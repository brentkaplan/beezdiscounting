# Converts answers from 5.5 trial probability discounting from Qualtrics template

Converts answers from 5.5 trial probability discounting from Qualtrics
template

## Usage

``` r
ans_pd(df)
```

## Arguments

- df:

  A dataframe containing all the columns

## Value

A dataframe with the ResponseId, index, and response (sc or lu).

## Examples

``` r
ans_pd(five.fivetrial_pd)
#> # A tibble: 21 × 3
#>    ResponseId index response
#>         <int> <chr> <chr>   
#>  1          1 I8    lu      
#>  2          1 I9    lu      
#>  3          1 I10   sc      
#>  4          1 I12   sc      
#>  5          1 I16   sc      
#>  6          2 I1    sc      
#>  7          2 I2    sc      
#>  8          2 I4    sc      
#>  9          2 I8    sc      
#> 10          2 I16   sc      
#> # ℹ 11 more rows
```
