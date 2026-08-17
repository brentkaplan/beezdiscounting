# Prepare data for the TMB mixed-effects discounting model

Operates on the FULL model frame from `.dd_validate_ip()` (canonical
`id`/`x`/`y` PLUS any retained factor/covariate columns). It
complete-cases that frame ONCE over all modeling columns (`id`, `x`,
`y`, and `extra_cols`), then derives the parallel arrays
`y`/`x`/`subject_id` from the same filtered rows so the later design
matrix, TMB arrays, and `fit$data` all share one row order
(row-coherence; mirrors beezdemand's `data_for_design`). Builds a
0-indexed `subject_id` aligned to `subject_levels` (the C++ template
indexes `u(subject_id, 0)` from 0).

## Usage

``` r
.dd_tmb_prepare_data(
  data,
  y_var = "y",
  x_var = "x",
  id_var = "id",
  extra_cols = NULL
)
```

## Arguments

- data:

  Data frame already coerced/clamped by `.dd_validate_ip()`, retaining
  factor/covariate columns alongside id/x/y.

- y_var, x_var, id_var:

  Character column names.

- extra_cols:

  Character vector of factor/covariate column names that must also be
  free of NA for a row to be kept (the union of `factors` and
  `continuous_covariates`). Defaults to none.

## Value

A list with `y`, `x`, `subject_id` (0-indexed integer),
`subject_levels`, `n_subjects`, `n_obs`, and `data` (the SINGLE
complete- cased model frame: canonical `id`, `x`, `y` first, then the
retained factor/covariate columns, in `data`'s row order).
