# Fitted values for a beezdiscounting_choice fit

Fitted values for a beezdiscounting_choice fit

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
fitted(object, level = c("subject", "population"), ...)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

- level:

  `"subject"` (default) or `"population"`.

- ...:

  Unused.

## Value

Numeric vector of fitted choice probabilities, length `nobs(object)`.
