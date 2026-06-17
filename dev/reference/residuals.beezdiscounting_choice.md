# Residuals for a beezdiscounting_choice fit

Residuals for a beezdiscounting_choice fit

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
residuals(
  object,
  type = c("response", "pearson"),
  level = c("subject", "population"),
  ...
)
```

## Arguments

- object:

  A `beezdiscounting_choice` object.

- type:

  `"response"` (default; raw `choice - p`) or `"pearson"` (divided by
  `sqrt(p * (1 - p))`).

- level:

  `"subject"` (default) or `"population"`.

- ...:

  Unused.

## Value

Numeric vector of residuals, length `nobs(object)`.
