# Converts answers from 5.5 trial delay discounting from Qualtrics template

Converts answers from 5.5 trial delay discounting from Qualtrics
template

## Usage

``` r
ans_dd(df)
```

## Arguments

- df:

  A dataframe containing all the columns

## Value

A dataframe with the ResponseId, index, and response (ss or ll).

## Examples

``` r
ans_dd(five.fivetrial_dd)
#> # A tibble: 22 × 3
#>    ResponseId index response
#>         <int> <chr> <chr>   
#>  1          1 I16   ll      
#>  2          1 I20   ll      
#>  3          1 I21   ss      
#>  4          1 I22   ss      
#>  5          1 I24   ss      
#>  6          2 I4    ll      
#>  7          2 I5    ss      
#>  8          2 I6    ss      
#>  9          2 I8    ss      
#> 10          2 I16   ss      
#> # ℹ 12 more rows
```
