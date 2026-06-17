# Run a single TMB optimization (nlminb or L-BFGS-B), family-agnostic

Ported verbatim from beezdemand `.tmb_run_optimizer` (renamed only).
Dispatches to `nlminb` (default) or `optim(method = "L-BFGS-B")` and
normalizes the return value so downstream code always sees identical
field names regardless of optimizer: `$par`, `$objective`,
`$convergence`, `$message` (guaranteed non-`NULL` character).

## Usage

``` r
.dd_tmb_run_optimizer(obj, start, tmb_control, user_specified, verbose)
```

## Arguments

- obj:

  TMB objective function object (from
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)), with
  `$fn`, `$gr`, and `$par`.

- start:

  Named numeric vector of starting parameter values.

- tmb_control:

  Merged control list; must have fields `optimizer`, `iter_max`,
  `eval_max`, `rel_tol`, `lower`, `upper`, `trace`.

- user_specified:

  Character vector of fields the user explicitly set in `tmb_control`
  (governs `trace` precedence).

- verbose:

  Integer verbosity level (0 = silent, 1 = progress, 2 = full trace).

## Value

A list with `opt` (list of `par`, `objective`, `convergence`, `message`)
and `warnings` (character vector).

## Details

Optimizer warnings are captured via `withCallingHandlers` and returned
in `$warnings`; errors yield `convergence = 99L, objective = Inf`.
