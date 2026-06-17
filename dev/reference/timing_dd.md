# Extract timing metrics from 5.5 trial delay discounting from Qualtrics template

Extract timing metrics from 5.5 trial delay discounting from Qualtrics
template

## Usage

``` r
timing_dd(df)
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
timing_dd(five.fivetrial_dd)
#> # A tibble: 22 × 7
#>    ResponseId index     q firstclick lastclick pagesubmit totalclicks
#>         <int> <chr> <dbl>      <dbl>     <dbl>      <dbl>       <dbl>
#>  1          1 I16       1      1.76      1.76        3.34           1
#>  2          1 I24       2      7.73      7.73        8.46           1
#>  3          1 I20       3      1.56      1.56        3.38           1
#>  4          1 I22       4      2.33      3.95        4.50           2
#>  5          1 I21       5      3.16      3.16        3.73           1
#>  6          2 I16       1      3.78      3.78        4.35           1
#>  7          2 I8        2      1.45      1.45        3.19           1
#>  8          2 I4        3      1.18      1.18        3.14           1
#>  9          2 I6        4      0.873     0.873       3.26           1
#> 10          2 I5        5      2.62      2.62        3.26           1
#> # ℹ 12 more rows
```
