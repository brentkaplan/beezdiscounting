# Reshape MCQ data long to wide

Reshape MCQ data long to wide

## Usage

``` r
long_to_wide_mcq(dat, q_col = "questionid", ans_col = "response")
```

## Arguments

- dat:

  Long format MCQ

- q_col:

  Name of the question column (default is "questionid")

- ans_col:

  Name of the answer column (default is "response")

## Value

Wide format data frame

## Examples

``` r
long_to_wide_mcq(mcq27)
#> # A tibble: 2 × 28
#>   subjectid   `1`   `2`   `3`   `4`   `5`   `6`   `7`   `8`   `9`  `10`  `11`
#>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1         1     0     0     0     1     1     0     1     1     0     0     1
#> 2         2     0     1     1     1     1     0     1     1     0     1     1
#> # ℹ 16 more variables: `12` <dbl>, `13` <dbl>, `14` <dbl>, `15` <dbl>,
#> #   `16` <dbl>, `17` <dbl>, `18` <dbl>, `19` <dbl>, `20` <dbl>, `21` <dbl>,
#> #   `22` <dbl>, `23` <dbl>, `24` <dbl>, `25` <dbl>, `26` <dbl>, `27` <dbl>
```
