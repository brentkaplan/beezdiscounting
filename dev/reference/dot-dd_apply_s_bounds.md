# Impose wide optimizer bounds on log_s for the 2-parameter equations

Mirrors
[`.dd_apply_phi_floor()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_apply_phi_floor.md):
adds `log_s` lower/upper entries to `tmb_control$lower`/`$upper` (wide
numerical guards `[log(0.05), log(20)]`) only when `has_s` and only when
the caller has not already set a `log_s` bound (user value wins). For
1-parameter equations `log_s` is mapped out of the optimizer vector, so
no bound is added (and `.expand_bounds` would otherwise warn about an
unknown parameter).

## Usage

``` r
.dd_apply_s_bounds(tmb_control, has_s)
```

## Arguments

- tmb_control:

  Merged control list (may carry user `lower`/`upper`).

- has_s:

  Logical; `TRUE` for green-myerson / rachlin.

## Value

The (possibly augmented) `tmb_control` list.
