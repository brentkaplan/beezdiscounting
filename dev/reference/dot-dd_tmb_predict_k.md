# Rebuild the per-row log-k linear predictor from newdata

Reconstructs the fixed-effect design from `newdata` using the stored RHS
formula and contrasts, multiplies by `beta_k`, and (for
`level = "subject"`) adds each subject's `sigma_u * u_i` deviate looked
up by id. Returns per-row `k` on the natural scale.

## Usage

``` r
.dd_tmb_predict_k(object, newdata, level = c("subject", "population"))
```

## Arguments

- object:

  A `beezdiscounting_tmb` fit.

- newdata:

  Data frame with the x var, factor/covariate cols, and (for
  `level = "subject"`) the id var.

- level:

  `"subject"` or `"population"`.

## Value

Numeric vector of per-row `k`, length `nrow(newdata)`.

## Details

Unseen factor levels in `newdata` are rejected with a clear error rather
than silently zero-padded (R2 fix).
