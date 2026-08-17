# Generate default starting values for the TMB discounting model

`beta_k` intercept is data-driven: invert the discounting function at
the median `y` observed at the minimum delay. `log_sigma_u = log(0.5)`.
The generic auxiliary scalar `log_aux` starts at `log(8)` for SLT-beta
(population precision phi) and `log(0.1)` for Gaussian (residual SD).

## Usage

``` r
.dd_tmb_default_starts(prepared, design, family, equation = "mazur")
```

## Arguments

- prepared:

  Output from
  [`.dd_tmb_prepare_data()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_prepare_data.md).

- design:

  Output from
  [`.dd_tmb_build_design()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_build_design.md).

- family:

  One of "sltb", "gaussian".

- equation:

  One of "mazur", "exponential", "green-myerson", "rachlin". The
  start-value inversion handles "mazur"/"exponential"; "green-myerson"/
  "rachlin" use the mazur start branch (`log_s` then starts at 0).

## Value

A parameters list: `beta_k` (length `ncol(X)`), `log_sigma_u`,
`log_aux`, `log_s` (always present; held fixed for 1-parameter
equations), `u` (matrix `n_subjects` x 1), plus the joint 2-RE blocks
`log_sd_re` (length 2), `cor_re` (length 1), and `b` (matrix
`n_subjects` x 2). The shared C++ template declares every block; the map
fixes the inactive ones per `n_re`.
