# Check for Unsystematic Data Violations

This function checks a dataset for violations of two criteria commonly
used to identify unsystematic delay-discounting data:

- Criterion 1: Any subsequent value of `y` exceeds the previous value by
  more than a specified proportion of the larger later reward (`ll`).

- Criterion 2: The last value of `y` is not at least a specified
  proportion less than the first value of `y`.

## Usage

``` r
check_unsystematic(dat, ll = 1, c1 = 0.2, c2 = 0.1)
```

## Arguments

- dat:

  A data frame containing the delay-discounting data. It must have at
  least two columns:

  - `id`: A unique identifier for the data set.

  - `y`: The indifference points to be analyzed.

  - `x` (optional): Delay values; when present, each subject's points
    are ordered by delay first.

- ll:

  A numeric value representing the larger later reward. Default is 1.

- c1:

  A numeric value for the threshold proportion for Criterion 1. Default
  is 0.2.

- c2:

  A numeric value for the threshold proportion for Criterion 2. Default
  is 0.1.

## Value

A tibble with the following columns:

- `id`: The unique identifier for the data set.

- `c1_pass`: Logical value indicating whether Criterion 1 was passed.

- `c2_pass`: Logical value indicating whether Criterion 2 was passed.

## Details

The criteria are applied separately to each `id`, so a data frame
containing several subjects returns one row per subject. When an `x`
(delay) column is present, each subject's points are ordered by delay
before the criteria are applied.

## Examples

``` r
data <- tibble::tibble(
  id = c(rep("P1", 6)),
  x = c(1, 7, 30, 90, 180, 365), # delays
  y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05) # indifference points
)
check_unsystematic(data, ll = 1, c1 = 0.2, c2 = 0.1)
#> # A tibble: 1 × 3
#>   id    c1_pass c2_pass
#>   <chr> <lgl>   <lgl>  
#> 1 P1    TRUE    TRUE   
```
