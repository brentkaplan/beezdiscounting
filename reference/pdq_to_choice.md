# Convert PDQ responses to a trial-level choice frame

Reshapes long-form Probability Discounting Questionnaire (PDQ) responses
into a per-trial certain-versus-risky choice frame, joining each
`questionid` to the canonical item design (certain amount, probabilistic
amount, probability, and odds against `theta = (1 - p)/p`) from Madden,
Petry, & Johnson (2009; see
[`get_lookup_table()`](https://brentkaplan.github.io/beezdiscounting/reference/get_lookup_table.md)).

## Usage

``` r
pdq_to_choice(
  responses,
  id_var = "subjectid",
  question_var = "questionid",
  response_var = "response"
)
```

## Arguments

- responses:

  Long-form data frame with one row per PDQ item per subject, holding
  the columns named by `id_var`, `question_var`, and `response_var`.
  `response` is `0` for the smaller guaranteed reward and `1` for the
  larger risky reward.

- id_var, question_var, response_var:

  Column names in `responses` for the subject id, PDQ question id
  (1-30), and the binary choice. Defaults match the bundled `pdq`
  dataset (`"subjectid"`, `"questionid"`, `"response"`).

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with
columns `id` (character), `sc_amount`, `lu_amount`, `prob`, `theta`
(odds against winning), and `choice` (`0`/`1`, `1` = chose the risky
reward), in the input row order.

## Details

Unknown or non-coercible question ids raise an error rather than
silently producing unmatched rows. Ragged input is allowed (subjects
need not have all items) and `NA` responses are preserved. For the
strict scorer see
[`score_pdq()`](https://brentkaplan.github.io/beezdiscounting/reference/score_pdq.md).

## See also

[`score_pdq()`](https://brentkaplan.github.io/beezdiscounting/reference/score_pdq.md),
[`get_lookup_table()`](https://brentkaplan.github.io/beezdiscounting/reference/get_lookup_table.md),
[`mcq_to_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/mcq_to_choice.md)

## Examples

``` r
ch <- pdq_to_choice(pdq)
head(ch)
#> # A tibble: 6 × 6
#>   id    sc_amount lu_amount  prob theta choice
#>   <chr>     <dbl>     <dbl> <dbl> <dbl>  <dbl>
#> 1 1            20        80  0.1   9         0
#> 2 1            20        80  0.13  6.69      0
#> 3 1            20        80  0.17  4.88      0
#> 4 1            20        80  0.2   4         0
#> 5 1            20        80  0.25  3         0
#> 6 1            20        80  0.33  2.03      1
```
