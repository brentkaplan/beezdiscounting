# Discounting function value (mean indifference proportion) for k and x

Mazur: `mu = 1 / (1 + k * x)`; exponential: `mu = exp(-k * x)`;
Green-Myerson: `mu = (1 + k * x)^(-s)`; Rachlin:
`mu = 1 / (1 + k * x^s)` with `x = 0 -> mu = 1`. Guards mu to
`[1e-6, 1 - 1e-6]` to match the C++ template bounds. The two 2-parameter
forms reduce to Mazur at `s = 1`.

## Usage

``` r
.dd_discount_mu(k, x, equation, s = 1)
```

## Arguments

- k:

  Numeric vector of discount rates.

- x:

  Numeric vector of delays (same length as `k`).

- equation:

  Character: `"mazur"`, `"exponential"`, `"green-myerson"`, or
  `"rachlin"`.

- s:

  Numeric nonlinearity exponent (Green-Myerson / Rachlin). Default `1`
  so 1-parameter callers are unchanged.

## Value

Numeric vector of mu values clamped to `[1e-6, 1-1e-6]`.
