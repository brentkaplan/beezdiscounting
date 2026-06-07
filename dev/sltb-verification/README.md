# SLT-beta density verification (pre-implementation provenance)

These scripts verified the **Scale-Location-Truncated Beta (SLT-beta)** discounting
density *before* it was ported into beezdiscounting's TMB mixed-effects tier. They are kept
as provenance; their property checks are formalized as package tests in
`tests/testthat/test-sltb-density.R`. (`^dev$` is `.Rbuildignore`d — not shipped.)

## What was verified

`scale_location_trun_beta_mazur` (reference) implements, with `a = μφ`, `b = (1−μ)φ`,
`s = 1.0000001`, `l = 1e-8`:

```
log f(y|μ,φ) = lgamma(a+b) − lgamma(a) − lgamma(b)
             + (a−1)·log(y/s+l) + (b−1)·log(1−(y/s+l))
             − log(s) − log(pbeta(1/s+l, a, b) − pbeta(l, a, b))
```

### `verify_sltb.R` — mathematical/numerical properties (all PASS)
- Normalizer: numeric kernel integral == analytic `Z = pbeta(1/s+l,a,b) − pbeta(l,a,b)`
  to 1e-14; normalized density integrates to 1 to 1e-13.
- **`Z` is load-bearing** (`Z(a=.05,b=1)=0.604`) — the `−log(Z)` term must stay.
- Finite log-density at `y=0` and `y=1` for all shapes (the SLT design feature).
- Scale/location limit: `f_SLT == dbeta/Z` to 1e-5.
- Moments: MC `E ≈ s(μ−l)`, `Var ≈ s²μ(1−μ)/(φ+1)` to truncation order.
- Our log-density extraction == the reference function body to 1e-12.
- MLE recovery (log k, log φ, multi-start): Mazur cor=0.986, exponential cor=0.953.

### `empirical_brent.R` — real-data validation on Jarvis et al. (2019) (`Data_brent.csv`)
- N=126 records, 7 delays (1wk–8yr); IP grid includes exact 0 and 1.
- **28 subjects have boundary IPs; standard beta is non-finite on exactly those 28; SLT is
  finite on all 126 and converges 126/126.** This is the core case for SLT.
- SLT-k vs NLS-k: cor(log k)=0.96; median k≈0.007; geometric-mean k=0.00733
  (the tie-out target for `exp(β̂₀_pop)`).
- Data note: `Data_brent.csv` `PID` is non-unique (78/126) — use the row as the subject
  unit; a mixed model must not group duplicated PIDs as one subject.

## Audit findings folded into the design (`docs/specs/2026-06-06-...md`)
1. Reference `ll.temp * ifelse(IP∈[0,1],1,0)` mask is buggy (`0*NaN=NaN`) → **filter rows, never mask.**
2. Reference optim is unconstrained → **estimate `log k`, `log φ`.**
3. Exponential μ underflows at long delays → **guard μ∈[1e-6,1−1e-6].**
4. Fixed start fails for exponential → **multi-start is essential.**
5. Reference `var.y` drops `s²` and truncation → derive from SLT moments / document.
6. Constants: paper prose says `10^−8.5`/`10^−9`; **code uses `s=1.0000001`, `l=1e-8`** (use code).

## References
- Kim, Koffarnus & Franck (2024). *Thinking Inside the Bounds...* (PMC11294315).
- Kim, Kaplan, Koffarnus & Franck (2025). *Scale-Location-Truncated Beta Regression.* arXiv:2509.13167.
- Reference code (collaborators', not vendored here): the `paper-code-kmg-slt_beta` repo;
  `R_code/beta functions.R`.
