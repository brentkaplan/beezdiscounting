# Calculate proportion of guaranteed choices at each PDQ h rank

The PDQ analog of
[`prop_ss()`](https://brentkaplan.github.io/beezdiscounting/reference/prop_ss.md):
pools all subjects' responses and reports, for each of the 10 h ranks
(one item per block per rank), the proportion choosing the smaller
guaranteed reward. `1 - prop_sc` is the risky-choice proportion at that
rank.

## Usage

``` r
prop_sc(dat)
```

## Arguments

- dat:

  Dataframe (longform) with subjectid, questionid (1-30), and response
  (0 for the guaranteed reward and 1 for the risky reward)

## Value

Dataframe with proportion of guaranteed choices at each h rank

## Examples

``` r
prop_sc(pdq)
#> # A tibble: 10 × 2
#>    h_rank prop_sc
#>     <int>   <dbl>
#>  1      1     0.5
#>  2      2     0.5
#>  3      3     0.5
#>  4      4     0.5
#>  5      5     0.5
#>  6      6     0  
#>  7      7     0  
#>  8      8     0  
#>  9      9     0  
#> 10     10     0  
```
