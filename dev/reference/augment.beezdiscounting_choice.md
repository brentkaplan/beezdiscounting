# Augment a beezdiscounting_choice model

Returns the data used for fitting (or `newdata`) as a tibble with three
diagnostic columns appended:

- `.fitted`: subject-conditional fitted choice probability P(LL).

- `.resid`: raw (response) residual `choice - .fitted`.

- `.std_resid`: Pearson residual `.resid / sqrt(p * (1 - p))`.

## Usage

``` r
# S3 method for class 'beezdiscounting_choice'
augment(x, newdata = NULL, ...)
```

## Arguments

- x:

  A `beezdiscounting_choice` object.

- newdata:

  Optional data frame. `NULL` uses the fitting data.

- ...:

  Unused.

## Value

A tibble with the same rows as the data plus `.fitted`, `.resid`, and
`.std_resid`.
