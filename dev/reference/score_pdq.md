# Score the 30-item Probability Discounting Questionnaire (PDQ)

Scores the Madden, Petry, & Johnson (2009) PDQ: three 10-item blocks
(block 1: \$20 for sure vs. a chance of \$80; block 2: \$40 vs. \$100;
block 3: \$40 vs. \$60), each an ascending ladder of h values at
indifference under the hyperbolic odds model `V = A / (1 + h * theta)`,
`theta = (1 - p) / p` (Rachlin, Raineri, & Cross, 1991). Each block is
scored independently by the same consistency-maximization algorithm as
the MCQ scorers; this implementation reproduces the Gray et al. (2016)
scoring syntax lookup tables for every possible response pattern.

## Usage

``` r
score_pdq(
  dat = dat,
  impute_method = "none",
  round = 6,
  random = FALSE,
  trans = "none",
  return_data = FALSE,
  verbose = FALSE
)
```

## Arguments

- dat:

  Dataframe (longform) with subjectid, questionid (1-30), and response
  (0 for the smaller guaranteed reward, 1 for the larger risky reward)

- impute_method:

  One of: "none", "ggm", "GGM", "inn", "INN"

- round:

  Numeric specifying number of decimal places (passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html))

- random:

  Boolean whether to insert a random draw (0 or 1) for NAs. Default is
  FALSE

- trans:

  Transformation to apply to h values: "none", "log", or "ln". Default
  is "none"

- return_data:

  Boolean whether to return the original data and new imputed responses.
  Default is FALSE.

- verbose:

  Boolean whether to print subject and question ids pertaining to
  missing data. Default is FALSE.

## Value

If `return_data = FALSE` (default), a summary data frame with one row
per subject: pooled `overall_h` (see Details), per-block h (`block1_h`,
`block2_h`, `block3_h`), their arithmetic mean (`mean_h`; Gray et al.'s
recommended composite) and geometric mean (`geomean_h`), pooled and
per-block consistency plus their mean (`composite_consistency`), and
pooled plus per-block proportions of risky choices (`block*_proportion`
is Gray et al.'s risky choice ratio). If `return_data = TRUE`, a list
with `results` and `data` (the input plus a `newresponse` column
reflecting any imputation).

## Details

Each subject's data must satisfy a strict contract: exactly one row per
canonical question id (30 of them; no duplicates, no unknown ids, none
missing) and responses coded 0, 1, or `NA` (numeric, logical, or
character/factor values that coerce to 0/1). Malformed input errors
rather than silently mis-scoring. Contrast with
[`pdq_to_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/pdq_to_choice.md)'s
lenient, ragged contract.

The published scoring (Madden et al., 2009; Gray et al., 2016) has no
overall 30-item ladder; blocks are scored separately, and `mean_h` is
Gray et al.'s recommended composite. `overall_h` and
`overall_consistency` are a beezdiscounting extension: all 30 items are
pooled into a single ascending ladder (exact-rational h order, ties
broken by question id; four item pairs tie exactly) and scored by the
same consistency-maximization algorithm with the repeat-last edge.

INN imputation groups items sharing an h rank (one item per block).
Ladders with remaining `NA` responses score `NA`; Gray et al. recommend
excluding subjects below 80% consistency on any block.

## References

Madden, G. J., Petry, N. M., & Johnson, P. S. (2009). Pathological
gamblers discount probabilistic rewards less steeply than matched
controls. *Experimental and Clinical Psychopharmacology, 17*(5),
283-290. [doi:10.1037/a0016806](https://doi.org/10.1037/a0016806)

Gray, J. C., Amlung, M. T., Palmer, A. A., & MacKillop, J. (2016).
Syntax for calculation of discounting indices from the monetary choice
questionnaire and probability discounting questionnaire. *Journal of the
Experimental Analysis of Behavior, 106*(2), 156-163.
[doi:10.1002/jeab.221](https://doi.org/10.1002/jeab.221)

## Examples

``` r
score_pdq(pdq)
#>   subjectid overall_h block1_h block2_h block3_h   mean_h geomean_h
#> 1         1  1.224745 1.215571 1.224745 1.233988 1.224768  1.224745
#> 2         2  0.329268 0.333333 0.329268 0.333333 0.331978  0.331973
#>   overall_consistency block1_consistency block2_consistency block3_consistency
#> 1                   1                  1                  1                  1
#> 2                   1                  1                  1                  1
#>   composite_consistency overall_proportion block1_proportion block2_proportion
#> 1                     1                0.5               0.5               0.5
#> 2                     1                1.0               1.0               1.0
#>   block3_proportion impute_method
#> 1               0.5          none
#> 2               1.0          none
```
