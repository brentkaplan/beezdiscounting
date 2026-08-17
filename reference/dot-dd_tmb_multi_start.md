# Multi-start optimization for the TMB discounting model

Builds 3 starting sets (data-driven, low-k, high-k) and runs each
through
[`.dd_tmb_run_optimizer()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_run_optimizer.md).
Keeps the lowest **finite** nll among candidates passing a NaN/blowup
sanity check on `beta_k[1]`.

## Usage

``` r
.dd_tmb_multi_start(
  tmb_data,
  start_values,
  tmb_control,
  user_specified,
  verbose,
  has_s = FALSE,
  n_re = 1L,
  covariance = "pdSymm"
)
```

## Arguments

- tmb_data:

  TMB data list from
  [`.dd_tmb_build_tmb_data()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_build_tmb_data.md).

- start_values:

  Default starting list from
  [`.dd_tmb_default_starts()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_default_starts.md).

- tmb_control:

  Merged control list (may carry user `lower`/`upper`).

- user_specified:

  Character vector of user-set `tmb_control` fields.

- verbose:

  Integer verbosity level.

- has_s:

  Logical; `TRUE` for the 2-parameter equations.

- n_re:

  Number of random-effect intercepts (`1L` or `2L`). Defaults to `1L` so
  existing positional callers stay valid.

- covariance:

  `"pdSymm"` or `"pdDiag"` (2-RE only).

## Value

A list with elements `obj`, `opt`, `nll`, `start_idx`, `opt_warnings`.

## Details

**Degenerate-optimum guard (phi-\>0 / k-\>inf, spec section 4.8):** For
`family == "sltb"` (i.e. `tmb_data$family == 0L`) a **lower bound on
`log_aux`** of `log(.dd_phi_min)` (= `log(0.1)`) is inserted into
`tmb_control$lower` before every optimizer call so the optimizer can
never walk into the degenerate phi-\>0 optimum. The bound is
overridable: if the caller already provides a `log_aux` entry in
`tmb_control$lower` that value wins (merged before the floor check).
Genuine low-precision fits (true phi above the floor, e.g. phi = 2) are
retained — the bound only blocks the pathological phi-\>0 sink. There is
**no post-hoc phi-based rejection**.

If every candidate trips the blowup guard the lowest-nll fit overall is
returned with a warning. If every candidate fails entirely (error or
non-finite nll) the function stops.
