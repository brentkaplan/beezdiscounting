# Score one subject's 27-item MCQ

Score one subject's 27-item MCQ

## Usage

``` r
score_one_mcq27(dat, impute_method = "none", round = 6)
```

## Arguments

- dat:

  One subject's 27 items from the MCQ

- impute_method:

  One of: "none", "ggm", "GGM", "inn", "INN"

- round:

  Numeric specifying number of decimal places (passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html))

## Value

Vector with scored 27-item MCQ metrics

## Examples

``` r
beezdiscounting:::score_one_mcq27(mcq27[mcq27$subjectid %in% 1, ])
#>             overall_k               small_k              medium_k 
#>              0.065212              0.025565              0.063690 
#>               large_k             geomean_k   overall_consistency 
#>              0.064947              0.047289              0.962963 
#>     small_consistency    medium_consistency     large_consistency 
#>              1.000000              1.000000              1.000000 
#> composite_consistency    overall_proportion      small_proportion 
#>              1.000000              0.259259              0.333333 
#>     medium_proportion      large_proportion 
#>              0.222222              0.222222 
```
