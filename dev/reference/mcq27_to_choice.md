# Convert 27-item MCQ responses to a trial-level choice frame

Reshapes long-form 27-item Monetary Choice Questionnaire (MCQ) responses
into the per-trial smaller-sooner versus larger-later choice frame
consumed by
[`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md),
joining each `questionid` to the canonical Kirby, Petry, & Bickel (1999)
item design (immediate amount, delayed amount, delay) bundled in the
lookup table (see
[`get_lookup_table()`](https://brentkaplan.github.io/beezdiscounting/reference/get_lookup_table.md)).

## Usage

``` r
mcq27_to_choice(
  responses,
  id_var = "subjectid",
  question_var = "questionid",
  response_var = "response"
)
```

## Arguments

- responses:

  Long-form data frame with one row per MCQ item per subject, holding
  the columns named by `id_var`, `question_var`, and `response_var`.
  `response` is `0` for the smaller-immediate reward (SIR/SS) and `1`
  for the larger-delayed reward (LDR/LL) – the same coding
  [`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md)
  expects, so no recoding is applied.

- id_var, question_var, response_var:

  Column names in `responses` for the subject id, MCQ question id
  (1-27), and the binary choice. Defaults match the bundled `mcq27`
  dataset (`"subjectid"`, `"questionid"`, `"response"`).

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with
columns `id` (character), `ss_amount`, `ll_amount`, `delay` (days), and
`choice` (`0`/`1`, `1` = chose LL), in the input row order. Ready to
pass to
[`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md).

## Details

Unknown or non-coercible question ids raise an error rather than
silently producing unmatched rows. Ragged input is allowed – subjects
need not have all 27 items – and `NA` responses are preserved (they are
complete-cased by
[`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md)).
For the strict 27-item scorer see
[`score_mcq27()`](https://brentkaplan.github.io/beezdiscounting/reference/score_mcq27.md).

## See also

[`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md),
[`score_mcq27()`](https://brentkaplan.github.io/beezdiscounting/reference/score_mcq27.md),
[`get_lookup_table()`](https://brentkaplan.github.io/beezdiscounting/reference/get_lookup_table.md)

## Examples

``` r
ch <- mcq27_to_choice(mcq27)
head(ch)
#> # A tibble: 6 × 5
#>   id    ss_amount ll_amount delay choice
#>   <chr>     <dbl>     <dbl> <dbl>  <dbl>
#> 1 1            54        55   117      0
#> 2 1            55        75    61      0
#> 3 1            19        25    53      0
#> 4 1            31        85     7      1
#> 5 1            14        25    19      1
#> 6 1            47        50   160      0
# feeds directly into the structural choice model (requires TMB):
# fit_dd_choice(ch, mode = "structural", equation = "mazur")
```
