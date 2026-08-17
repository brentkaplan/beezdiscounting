# Extract timing metrics from 5.5 trial probability discounting from Qualtrics template

Extract timing metrics from 5.5 trial probability discounting from
Qualtrics template

## Usage

``` r
timing_pd(df)
```

## Arguments

- df:

  A dataframe containing all the columns

## Value

A dataframe with ResponseId, indexes, values and timing

## Details

Currently assumes the attending questions are present and labeled
"Attend-LL" and "Attend-SS"

## Examples

``` r
timing_pd(five.fivetrial_pd)
#> # A tibble: 21 × 7
#>    ResponseId index     q firstclick lastclick pagesubmit totalclicks
#>         <int> <chr> <dbl>      <dbl>     <dbl>      <dbl>       <dbl>
#>  1          1 I16       1       3.98      3.98       5.18           1
#>  2          1 I8        2       4.01      4.01       4.76           1
#>  3          1 I12       3       2.06      2.06       3.25           1
#>  4          1 I10       4       1.52      1.52       3.02           1
#>  5          1 I9        5       2.25      2.95       3.74           2
#>  6          2 I16       1       2.87      2.87       3.88           1
#>  7          2 I8        2       3.74      3.74       4.86           1
#>  8          2 I4        3       1.16      1.16       6.36           1
#>  9          2 I2        4       3.06      3.06       5.41           1
#> 10          2 I1        5       2.05      2.05       5.10           1
#> # ℹ 11 more rows
```
