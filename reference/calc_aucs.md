# Calculate Area-Under-the-Curve (AUC) Metrics for Delay Discounting Data

Calculates three area-under-the-curve (AUC) metrics for
delay-discounting data: regular AUC (using raw delays), log10 AUC (using
logarithmically scaled delays), and ordinal AUC (using ordinally scaled
delays). The three differ only in how the delays are scaled before the
area is computed. Metrics are computed separately for each `id`, so a
data frame with several subjects returns one row per subject.

## Usage

``` r
calc_aucs(dat)
```

## Arguments

- dat:

  A data frame containing delay discounting data. It must include the
  following columns:

  - `id`: Participant or group identifier.

  - `x`: Delay values (e.g., in days).

  - `y`: Indifference points as a proportion of the larger later reward
    (the area is normalized by the delay range times 1, so `y` must be
    on the 0–1 scale).

  Missing `y` or `x` values and duplicate delays within a subject are
  errors.

## Value

A tibble with the following columns:

- `id`: The participant or group identifier.

- `auc_regular`: The regular AUC, calculated using the raw delay values.

- `auc_log10`: The log10 AUC, calculated on `log10(x + 1)` delays scaled
  to their maximum.

- `auc_ord`: The ordinal AUC, calculated using ordinally scaled delay
  values.

## Details

Each area is the trapezoidal area under the indifference points divided
by the width of the delay axis, so a subject who does not discount
scores 1. `auc_log10` transforms delays as `log10(x + 1)` (the `+ 1`
keeps a zero delay finite, the convention usually attributed to Borges
et al., 2016) and rescales them by their maximum. Because of the `+ 1`,
`auc_log10` depends on the delay unit: the same series scores
differently in days and in weeks, so compare `auc_log10` only across
data recorded in one delay unit. `auc_regular` and `auc_ord` are
unit-free.

## References

Borges, A. M., Kuang, J., Milhorn, H., & Yi, R. (2016). An alternative
approach to calculating area-under-the-curve (AUC) in delay discounting
research. *Journal of the Experimental Analysis of Behavior, 106*,
145–155. [doi:10.1002/jeab.219](https://doi.org/10.1002/jeab.219)

## Examples

``` r
# Example data
data <- data.frame(
  id = rep("P1", 6),
  x = c(1, 7, 30, 90, 180, 365),
  y = c(0.8, 0.5, 0.3, 0.2, 0.1, 0.05)
)

# Calculate AUC metrics for a single participant
calc_aucs(data)
#> # A tibble: 1 × 4
#>   id    auc_regular auc_log10 auc_ord
#>   <chr>       <dbl>     <dbl>   <dbl>
#> 1 P1          0.152     0.359   0.305
```
