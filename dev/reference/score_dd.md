# Score 5.5 trial delay discounting from Qualtrics template

Score 5.5 trial delay discounting from Qualtrics template

## Usage

``` r
score_dd(df)
```

## Arguments

- df:

  A dataframe containing all the columns

## Value

A dataframe with id, indexes, response, k value, and effective delay 50.

## Details

Currently assumes the attending questions are present and labeled
"Attend-LL" and "Attend-SS"

## Examples

``` r
score_dd(five.fivetrial_dd)
#> # A tibble: 4 × 6
#>   ResponseId index    response     kval attentionflag    ed50
#>        <int> <chr>    <chr>       <dbl> <chr>           <dbl>
#> 1          1 I21      ss        0.00671 No            149.   
#> 2          2 I5       ss        4.90    No              0.204
#> 3          3 AttendSS ss       NA       Yes            NA    
#> 4          4 AttendLL ll       NA       Yes            NA    
```
