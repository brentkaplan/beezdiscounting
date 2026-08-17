# Score one ascending ladder by consistency maximization

The Kirby-style scoring core shared by the MCQ scorers (overall +
magnitude ladders) and the PDQ block scorer. Extracted verbatim from
[`score_one_mcq()`](https://brentkaplan.github.io/beezdiscounting/reference/score_one_mcq.md)'s
magnitude loop; reproduces Gray et al.'s (2016) published PDQ lookup
tables for all 1024 response patterns per block.

## Usage

``` r
.score_ladder(resp, vals, edge)
```

## Arguments

- resp:

  0/1/NA responses in ladder order (1 = chose the larger delayed/risky
  option). Any `NA` makes all three outputs `NA` – impute upstream.

- vals:

  Ascending indifference values (k or h), one per item; must be the same
  length as `resp`.

- edge:

  Value appended past the steepest item (the overall-ladder edge
  constant, or `vals[length(vals)]` for the repeat-last convention).

## Value

`list(value, consistency, proportion)`.
