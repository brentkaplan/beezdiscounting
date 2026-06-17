# Build the TMB `map` for the discounting shape + random-effect blocks

Generalizes the old `has_s`-only map to cover the joint 2-RE blocks. The
shared C++ template declares EVERY parameter (`log_s`,
`log_sigma_u`/`u`, `log_sd_re`/`cor_re`/`b`); this map fixes the
inactive ones so a given fit only estimates the blocks it should:

## Usage

``` r
.dd_tmb_build_map(has_s, n_re = 1L, covariance = "pdSymm")
```

## Arguments

- has_s:

  Logical; `TRUE` for green-myerson / rachlin.

- n_re:

  Number of random-effect intercepts (`1L` or `2L`).

- covariance:

  `"pdSymm"` (cor_re free) or `"pdDiag"` (cor_re fixed); only consulted
  when `n_re == 2L`.

## Value

A (possibly empty) named list of `factor(NA)` map entries.

## Details

- `log_s` is held fixed (`factor(NA)`) for the 1-parameter equations
  (`!has_s`) and freed for the 2-parameter equations (`has_s`).

- `n_re == 1L` (single intercept on `log k`): fix `log_sd_re`, `cor_re`
  (and `b` is fixed by
  [`.dd_tmb_finalize_map()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_finalize_map.md),
  which needs the starts to size the factor).

- `n_re == 2L` (joint `(log k, log phi)`): fix `log_sigma_u` (and `u`
  via
  [`.dd_tmb_finalize_map()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_finalize_map.md));
  fix `cor_re` iff `covariance == "pdDiag"`.

The returned value is finalized by
[`.dd_tmb_finalize_map()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_finalize_map.md)
and passed verbatim to EVERY
[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html) call so
a fit can never estimate an unidentified block (singular Hessian, wrong
df).
