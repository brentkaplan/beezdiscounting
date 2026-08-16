# Internal instrument registry

Serves the per-instrument design table and ladder metadata for the
questionnaire scorers: the 27-item MCQ (Kirby, Petry, & Bickel, 1999),
the 21-item MCQ (Kirby & Marakovic, 1996), and the 30-item PDQ (Madden,
Petry, & Johnson, 2009).

## Usage

``` r
.instrument_registry(instrument)
```

## Arguments

- instrument:

  One of "mcq27", "mcq21", "pdq".

## Value

List with `key`, `items`, `table` (design table in ladder order),
`n_mag` (items per magnitude/block), `edge_k` (overall-ladder edge
constant, or `NA` when the overall ladder uses the repeat-last
convention instead), `overall_order_col` (name of the table column
giving each item's position in the pooled overall ladder, or `NULL` when
the table's row order already is that ladder), `value_col`, `rank_col`,
`param`, and `rank_labels`.

## Details

Edge conventions: MCQ overall ladders append 0.25 (27-item) / 0.1333
(21-item) past the steepest item (Kaplan et al., 2014 Excel scorers);
MCQ magnitude ladders and all PDQ block ladders repeat their last
indifference value (verified against all 3 x 1024 response patterns in
Gray et al.'s 2016 PDQ lookup tables). The PDQ's pooled overall ladder
(`overall_order_col = "overall_rank"`, repeat-last edge) is a
beezdiscounting extension – the published PDQ scoring (Madden et al.,
2009; Gray et al., 2016) scores blocks only.
