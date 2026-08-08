# Score one subject's MCQ (27- or 21-item)

Score one subject's MCQ (27- or 21-item)

## Usage

``` r
score_one_mcq(dat, reg, impute_method = "none", round = 6)
```

## Arguments

- dat:

  One subject's items from the MCQ

- reg:

  Registry list from
  [`.mcq_registry()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-mcq_registry.md)

- impute_method:

  One of: "none", "ggm", "GGM", "inn", "INN"

- round:

  Numeric specifying number of decimal places (passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html))

## Value

Vector with scored MCQ metrics
