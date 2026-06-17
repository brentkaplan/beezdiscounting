# Expected free-parameter blocks for a given map configuration

The exact set of optimizer parameter blocks left free by the map for a
`(has_s, n_re, covariance)` combination. Uses the raw optimizer name
`log_aux` (the rename to `log_phi`/`log_sigma_e` happens in
[`.dd_tmb_extract_estimates()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_extract_estimates.md)).

## Usage

``` r
.dd_expected_free_params(has_s, n_re, covariance)
```

## Arguments

- has_s:

  Logical; `TRUE` for the 2-parameter equations.

- n_re:

  Number of random-effect intercepts (`1L` or `2L`).

- covariance:

  `"pdSymm"` or `"pdDiag"` (2-RE only).

## Value

Character vector of expected free-block names.
