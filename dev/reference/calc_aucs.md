# Calculate Area-Under-the-Curve (AUC) Metrics for Delay Discounting Data

This function calculates three types of Area-Under-the-Curve (AUC)
metrics for delay discounting data: regular AUC (using raw delays),
log10 AUC (using logarithmically scaled delays), and ordinal AUC (using
ordinally scaled delays). These metrics provide different perspectives
on the rate of delay discounting. Metrics are computed separately for
each `id`, so a data frame with several subjects returns one row per
subject.

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

  - `y`: Indifference point values (e.g., subjective value of the
    delayed reward).

## Value

A tibble with the following columns:

- `id`: The participant or group identifier.

- `auc_regular`: The regular AUC, calculated using the raw delay values.

- `auc_log10`: The log10 AUC, calculated using logarithmically
  transformed delay values.

- `auc_ord`: The ordinal AUC, calculated using ordinally scaled delay
  values.

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
