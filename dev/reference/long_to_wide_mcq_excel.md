# Reshape MCQ data from long to wide (as used in the 21- and 27-Item Monetary Choice Questionnaire Automated Scorer)

Reshape MCQ data from long to wide (as used in the 21- and 27-Item
Monetary Choice Questionnaire Automated Scorer)

## Usage

``` r
long_to_wide_mcq_excel(dat, subj_col = "subjectid", ans_col = "response")
```

## Arguments

- dat:

  Long format MCQ data

- subj_col:

  Character column name of subject ids

- ans_col:

  Character column name of responses

## Value

Wide format MCQ data that can be used in the Excel Automated Scorers

## Examples

``` r
long_to_wide_mcq_excel(generate_data_mcq(2))
#> # A tibble: 27 × 3
#>    questionid   `1`   `2`
#>         <int> <int> <int>
#>  1          1     1     1
#>  2          2     1     1
#>  3          3     1     1
#>  4          4     1     0
#>  5          5     0     1
#>  6          6     1     0
#>  7          7     0     0
#>  8          8     0     0
#>  9          9     0     1
#> 10         10     1     0
#> # ℹ 17 more rows
```
