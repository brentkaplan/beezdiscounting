# Guard, clamp and floor activity at the fitted values

Counts, at the subject-level fitted values, the rows with a positive
delay whose raw mean lies outside the `[1e-6, 1 - 1e-6]` guard (`x = 0`
rows have mean 1 by definition and are excluded); for `k + s ~ 1`, the
subjects whose soft-clamped `s` differs materially from the latent value
(`.dd_s_clamp_active()`); for `k + phi ~ 1`, the subjects whose latent
`phi` is below the 0.1 floor. Reporting only (audit F-BZ4-2, F-BZ4-5).

## Usage

``` r
.dd_tmb_guard_info(object)
```

## Arguments

- object:

  A `beezdiscounting_tmb` fit (with `subject_pars`).

## Value

list(n_rows, mu_guard_lower, mu_guard_upper, and either
n_s_clamped_lower / n_s_clamped_upper or n_phi_floor when applicable).
