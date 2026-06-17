# Score 5.5 trial probability discounting from Qualtrics template

Score 5.5 trial probability discounting from Qualtrics template

## Usage

``` r
score_pd(df)
```

## Arguments

- df:

  A dataframe containing all the columns

## Value

A dataframe with id, indexes, response, h value, and effective
probability 50.

## Details

Currently assumes the attending questions are present and labeled
"Attend-LL" and "Attend-SS"

## Examples

``` r
score_pd(five.fivetrial_pd)
#> # A tibble: 4 × 7
#>   ResponseId index    response  hval attentionflag etheta50  ep50
#>        <int> <chr>    <chr>    <dbl> <chr>            <dbl> <dbl>
#> 1          1 I9       lu        7.44 No              0.134   88.1
#> 2          2 AttendSS lu       99    No              0.0101  99  
#> 3          3 I15      sc        1.60 No              0.624   61.6
#> 4          4 I9       lu        7.44 No              0.134   88.1
```
