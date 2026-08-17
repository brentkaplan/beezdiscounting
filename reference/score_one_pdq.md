# Score one subject's PDQ

Score one subject's PDQ

## Usage

``` r
score_one_pdq(dat, reg, impute_method = "none", round = 6)
```

## Arguments

- dat:

  One subject's 30 PDQ items in long form

- reg:

  Registry list from `.instrument_registry("pdq")`

- impute_method:

  One of: "none", "ggm", "GGM", "inn", "INN"

- round:

  Numeric specifying number of decimal places (passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html))

## Value

Named vector with scored PDQ metrics
