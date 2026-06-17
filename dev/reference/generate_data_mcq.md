# Generate fake MCQ data

Generate fake MCQ data

## Usage

``` r
generate_data_mcq(n_ids = 100, n_items = 27, seed = 1234, prop_na = 0)
```

## Arguments

- n_ids:

  Number of subjectids

- n_items:

  Number of trials

- seed:

  Random seed

- prop_na:

  Proportion of NAs in the entire data set

## Value

Dataframe of subjectid, questionid, and response

## Examples

``` r
generate_data_mcq(n_ids = 2, n_items = 27, prop_na = .01)
#>    subjectid questionid response
#> 1          1          1        1
#> 2          1          2        1
#> 3          1          3        1
#> 4          1          4        1
#> 5          1          5        0
#> 6          1          6        1
#> 7          1          7        0
#> 8          1          8        0
#> 9          1          9        0
#> 10         1         10        1
#> 11         1         11        1
#> 12         1         12        1
#> 13         1         13        1
#> 14         1         14        0
#> 15         1         15        1
#> 16         1         16        1
#> 17         1         17        1
#> 18         1         18        0
#> 19         1         19        1
#> 20         1         20        1
#> 21         1         21        1
#> 22         1         22        1
#> 23         1         23        1
#> 24         1         24        1
#> 25         1         25        1
#> 26         1         26        1
#> 27         1         27        0
#> 28         2          1        1
#> 29         2          2        1
#> 30         2          3        1
#> 31         2          4        0
#> 32         2          5        1
#> 33         2          6        0
#> 34         2          7        0
#> 35         2          8        0
#> 36         2          9        1
#> 37         2         10        0
#> 38         2         11        1
#> 39         2         12        1
#> 40         2         13        0
#> 41         2         14        1
#> 42         2         15        0
#> 43         2         16       NA
#> 44         2         17        1
#> 45         2         18        1
#> 46         2         19        0
#> 47         2         20        0
#> 48         2         21        0
#> 49         2         22        0
#> 50         2         23        1
#> 51         2         24        1
#> 52         2         25        1
#> 53         2         26        1
#> 54         2         27        1
```
