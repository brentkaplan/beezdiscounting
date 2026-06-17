# Impose the SLT-beta phi floor as a log_aux lower bound

For `family == "sltb"` the degenerate phi-\>0 / k-\>inf optimum (spec
section 4.8) is blocked by a **lower bound on `log_aux`** of
`log(.dd_phi_min)` inserted into `tmb_control$lower`. This must apply to
BOTH optimizer paths (single fit and multi-start), so the logic lives in
one shared helper. The bound is overridable: if the caller already
supplies a `log_aux` entry in `tmb_control$lower`, their value wins (the
floor is not added).

## Usage

``` r
.dd_apply_phi_floor(tmb_control, family)
```

## Arguments

- tmb_control:

  Merged control list (may carry user `lower`/`upper`).

- family:

  One of "sltb", "gaussian".

## Value

The (possibly augmented) `tmb_control` list.
