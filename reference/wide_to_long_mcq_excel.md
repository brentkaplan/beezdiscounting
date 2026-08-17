# Reshape MCQ data from wide (as used in the 21- and 27-Item Monetary Choice Questionnaire Automated Scorer) to long

Reshape MCQ data from wide (as used in the 21- and 27-Item Monetary
Choice Questionnaire Automated Scorer) to long

## Usage

``` r
wide_to_long_mcq_excel(dat)
```

## Arguments

- dat:

  Wide format MCQ data as used in the Excel Automated Scorers

## Value

Long format data frame

## Examples

``` r
wide_to_long_mcq_excel(long_to_wide_mcq_excel(generate_data_mcq(2)))
#> # A tibble: 54 × 3
#>    subjectid questionid response
#>    <chr>          <int>    <int>
#>  1 1                  1        1
#>  2 1                  2        1
#>  3 1                  3        1
#>  4 1                  4        1
#>  5 1                  5        0
#>  6 1                  6        1
#>  7 1                  7        0
#>  8 1                  8        0
#>  9 1                  9        0
#> 10 1                 10        1
#> # ℹ 44 more rows
```
