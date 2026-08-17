# Compute subject-specific discounting parameters

Reconstructs each subject's `k_i = exp(eta_i)` where
`eta_i = X_i %*% beta_k + sigma_u * u_i`, exactly matching the C++
template's non-centered predictor
(`log_k_i = X.row(i)*beta_k + sigma_u * u(subj,0)`). `X_i` is subject
`i`'s fixed-effect design ROW. For between-subject factors/covariates
the design is constant within a subject, so the first design-matrix row
for each subject is used. This makes per-subject `k`
factor/covariate-correct: subjects in non-reference groups pick up their
group's `beta_k` contribution rather than the reference-group intercept.

## Usage

``` r
.dd_tmb_compute_subject_pars(
  coefficients,
  u_hat,
  subject_levels,
  design_X,
  subject_id,
  equation,
  family,
  n_re = 1L,
  log_aux_name = "log_phi",
  Sigma = NULL,
  re2_target = 0L
)
```

## Arguments

- coefficients:

  Named coefficient vector (with `beta_k`, `log_sigma_u` for 1-RE, and
  `log_phi` or `log_sigma_e`).

- u_hat:

  Matrix of per-subject random effects: `n_subjects x 1` standardized
  `u` deviates (1-RE) or `n_subjects x 2` standardized `b` deviates
  (2-RE).

- subject_levels:

  Character vector of subject ids (length n_subjects).

- design_X:

  Fixed-effect design matrix (rows aligned with the prepared arrays);
  columns correspond 1:1 with the `beta_k` entries.

- subject_id:

  Integer 0-indexed subject map (length `n_obs`), aligned row-for-row
  with `design_X`. Used to locate each subject's first design row.

- equation:

  One of "mazur", "exponential" (reserved; k is equation-free).

- family:

  One of "sltb", "gaussian".

- n_re:

  Number of random-effect intercepts (`1L` or `2L`). Defaults to `1L` so
  existing positional callers stay valid.

- log_aux_name:

  Name of the auxiliary (precision) coefficient in `coefficients`
  (`"log_phi"` for sltb). Used only for `n_re == 2`.

- Sigma:

  Fitted 2x2 RE covariance (from
  [`.dd_tmb_extract_estimates()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_extract_estimates.md));
  required for `n_re == 2`.

- re2_target:

  Integer flag: `0L` for a phi-target 2-RE fit, `1L` for an s-target
  fit. Ignored when `n_re == 1L`. Defaults to `0L`.

## Value

data.frame: `id, u_i, k` (1-RE), `id, re_k, re_phi, k, phi` (phi-target
2-RE), or `id, re_k, re_s, k, s` (s-target 2-RE).

## Details

For a 1-RE fit the auxiliary scalar `phi` is population-level, so it is
not a subject-level parameter (the returned frame is `id, u_i, k`). For
a phi-target 2-RE fit (`k + phi ~ 1`) each subject's natural-scale
`re = L * b_hat` is reconstructed (`L = chol(Sigma)` lower-triangular;
the `sdreport` "random" block holds the STANDARDIZED `b_hat`, not `re`),
giving `id, re_k, re_phi, k, phi` with a per-subject `phi_i` floor
matching the kernel clamp. For an s-target 2-RE fit (`k + s ~ 1`) the
same reconstruction gives `id, re_k, re_s, k, s` with
`s_i = .dd_soft_clamp_s_log(log_s + re_s)` (approximately
`exp(log_s + re_s)` in the interior, smoothly saturating toward
`(0.05, 20)`).

## Note

Subject-level `k` assumes between-subject predictors: the first design
row per subject defines that subject's fixed-effect contribution. A
within-subject-varying covariate would make a single per-subject `k`
ill-defined; in that case the first observed row is used.
