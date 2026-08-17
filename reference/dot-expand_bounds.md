# Expand partial optimizer bounds to the full parameter vector

Ported verbatim from beezdemand `.expand_bounds`. Given a named vector
of user-specified bounds (possibly partial, possibly `NULL`), expands it
to the full length of `par_names` by filling unspecified positions with
`default_val`. Repeated parameter names (e.g. `beta_k`, `beta_k`) are
each filled when the corresponding name appears in `bounds`.

## Usage

``` r
.expand_bounds(bounds, par_names, default_val)
```

## Arguments

- bounds:

  Named numeric vector of user-specified bounds, or `NULL`.

- par_names:

  Character vector of optimizer parameter names (from `names(obj$par)`).
  May contain repeated names for vector parameters.

- default_val:

  Default bound value: `-Inf` for lower, `Inf` for upper.

## Value

Numeric vector of length `length(par_names)`.
