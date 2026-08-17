# Build the TMB data list for MixedDiscounting

Build the TMB data list for MixedDiscounting

## Usage

``` r
.dd_tmb_build_tmb_data(
  prepared,
  design,
  equation,
  family,
  n_re = 1L,
  re2_target = 0L
)
```

## Arguments

- prepared:

  Output from
  [`.dd_tmb_prepare_data()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_prepare_data.md).

- design:

  Output from
  [`.dd_tmb_build_design()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_build_design.md).

- equation:

  One of "mazur", "exponential", "green-myerson", "rachlin" (sets the
  kernel `eqn_type`).

- family:

  One of "sltb", "gaussian".

- n_re:

  Number of random-effect intercepts: `1L` (on `log k`) or `2L` (joint
  `(log k, second target)`). Defaults to `1L`.

- re2_target:

  Second-RE target for `n_re == 2L`: `0L` = the precision `log phi`
  (default), `1L` = the GM/Rachlin curvature `log s`. Ignored for
  `n_re == 1L`.

## Value

A list whose names match the C++ `DATA_*` macros: `model`, `y`, `x`,
`subject_id` (0-indexed integer), `X`, `eqn_type`, `family`, `n_obs`,
`n_subjects`, `n_re`, `re2_target`.
