# Fix the inactive random-effect matrix in the TMB map

`u` (1-RE) and `b` (2-RE) are alternately the integrated random block vs
a fixed-at-zero block. This appends the `factor(NA)` fix for whichever
matrix is inactive, sized from its starts. Applied right before each
[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

## Usage

``` r
.dd_tmb_finalize_map(map, starts, n_re)
```

## Arguments

- map:

  The map from
  [`.dd_tmb_build_map()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_build_map.md).

- starts:

  The starting-value list (for the matrix lengths).

- n_re:

  Number of random-effect intercepts (`1L` or `2L`).

## Value

The augmented map.
