# Reshape MCQ data wide to long

Reshape MCQ data wide to long

## Usage

``` r
wide_to_long_mcq(dat, items = 27)
```

## Arguments

- dat:

  Wide format MCQ assuming subject id is in column 1 and one column per
  item after it

- items:

  Number of MCQ questions; must equal `ncol(dat) - 1`

## Value

Long format data frame

## Details

Items are numbered by column position within each wide row, so the
result has `items` rows per input row whatever the subject ids are.
Repeated or missing subject ids are left for the scoring functions to
report.

## Examples

``` r
wide_to_long_mcq(long_to_wide_mcq(mcq27))
#> # A tibble: 54 × 3
#>    subjectid questionid response
#>        <dbl>      <int>    <dbl>
#>  1         1          1        0
#>  2         1          2        0
#>  3         1          3        0
#>  4         1          4        1
#>  5         1          5        1
#>  6         1          6        0
#>  7         1          7        1
#>  8         1          8        1
#>  9         1          9        0
#> 10         1         10        0
#> # ℹ 44 more rows
```
