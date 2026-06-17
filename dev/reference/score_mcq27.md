# Score 27-item MCQ

Score 27-item MCQ

## Usage

``` r
score_mcq27(
  dat = dat,
  impute_method = "none",
  round = 6,
  random = FALSE,
  trans = "none",
  return_data = FALSE,
  verbose = FALSE
)
```

## Arguments

- dat:

  Dataframe (longform) with subjectid, questionid, and response (0 for
  SIR/SS and 1 for LDR/LL)

- impute_method:

  One of: "none", "ggm", "GGM", "inn", "INN"

- round:

  Numeric specifying number of decimal places (passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html))

- random:

  Boolean whether to insert a random draw (0 or 1) for NAs. Default is
  FALSE

- trans:

  Transformation to apply to k values: "none", "log", or "ln". Default
  is "none"

- return_data:

  Boolean whether to return the original data and new imputed responses.
  Default is FALSE.

- verbose:

  Boolean whether to print subject and question ids pertaining to
  missing data. Default is FALSE.

## Value

Summary dataframe

## Examples

``` r
score_mcq27(mcq27)
#>   subjectid overall_k  small_k medium_k  large_k geomean_k overall_consistency
#> 1         1  0.065212 0.025565 0.063690 0.064947  0.047289            0.962963
#> 2         2  0.000399 0.000633 0.001589 0.000251  0.000632            0.962963
#>   small_consistency medium_consistency large_consistency composite_consistency
#> 1                 1                  1                 1                     1
#> 2                 1                  1                 1                     1
#>   overall_proportion small_proportion medium_proportion large_proportion
#> 1           0.259259         0.333333          0.222222         0.222222
#> 2           0.777778         0.777778          0.666667         0.888889
#>   impute_method
#> 1          none
#> 2          none
```
