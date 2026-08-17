# Assert the exact free-parameter block set at the map seam

Pins the free-block set for each `(n_re, covariance, has_s)` combination
so a leaked block (an unidentified parameter slipping into the
optimizer) or a missing block (a free coefficient accidentally mapped
out) is caught immediately rather than surfacing as a silent fit
pathology.

## Usage

``` r
.dd_assert_free_params(par_names, expected, context)
```

## Arguments

- par_names:

  Character vector of optimizer parameter names (`names(opt$par)`).

- expected:

  Character vector of expected free-block names.

- context:

  Short label for the error message.

## Value

Invisibly `TRUE`; stops on mismatch.
