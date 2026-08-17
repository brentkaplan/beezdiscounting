# Impose wide optimizer bounds on log_sd_re for a 2-RE fit

Adds lower/upper `log_sd_re` bounds (`[log(1e-3), log(5)]`) only when
`n_re == 2L` and the caller has not already set a `log_sd_re` bound.
Each bound is duplicated (the vector has length 2). Keeps Sigma
positive-definite and stops an absurd RE-variance estimate from
collapsing a boundary subject's `phi_i`.

## Usage

``` r
.dd_apply_re_bounds(tmb_control, n_re)
```

## Arguments

- tmb_control:

  Merged control list (may carry user `lower`/`upper`).

- n_re:

  Number of random-effect intercepts (`1L` or `2L`).

## Value

The (possibly augmented) `tmb_control` list.
