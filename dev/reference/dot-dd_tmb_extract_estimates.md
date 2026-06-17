# Extract estimates from a TMB discounting fit

Computes [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html),
gates on `isTRUE(sdr$pdHess)` (warn, never abort), and renames the
generic auxiliary scalar `log_aux` to `log_phi` (sltb) or `log_sigma_e`
(gaussian) in the returned `coefficients`/`se`.

## Usage

``` r
.dd_tmb_extract_estimates(
  obj,
  opt,
  n_subjects,
  family,
  has_s = FALSE,
  verbose = 1,
  n_re = 1L,
  re2_target = 0L
)
```

## Arguments

- obj:

  TMB objective object.

- opt:

  Normalized optimizer result (from
  [`.dd_tmb_run_optimizer()`](https://brentkaplan.github.io/beezdiscounting/reference/dot-dd_tmb_run_optimizer.md)).

- n_subjects:

  Integer number of subjects.

- family:

  One of "sltb", "gaussian".

- has_s:

  Logical; `TRUE` for the 2-parameter equations (green-myerson /
  rachlin), where `log_s` is a free coefficient.

- verbose:

  Integer verbosity.

- n_re:

  Number of random-effect intercepts (`1L` or `2L`). Defaults to `1L` so
  existing positional callers stay valid.

- re2_target:

  Second-RE target (`0L` = phi, `1L` = s); selects the fitted `Sigma`
  dimnames (`(k, phi)` vs `(k, s)`) for `n_re == 2L`. Defaults to `0L`.

## Value

list(coefficients, se, sdr, variance_components, u_hat, hessian_pd,
Sigma). `u_hat` is the per-subject random-effect block: a 1-column
matrix of standardized `u` deviates for `n_re == 1`, an `n_subjects x 2`
matrix of standardized `b` deviates for `n_re == 2`. `Sigma` is the
fitted 2x2 RE covariance for `n_re == 2` (NULL otherwise).
