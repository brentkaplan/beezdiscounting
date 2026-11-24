# Modeling delay-discounting indifference points with bounded error distributions

## Why indifference points need a bounded error distribution

A delay-discounting task asks a person to trade a smaller-sooner reward
against a larger-later one across a range of delays. The **indifference
point** (IP) at each delay is the smaller-sooner amount that feels
equivalent to the later reward, expressed as a *proportion of that later
reward*. By construction an IP lives on the closed interval \[0, 1\]:

- y = 1: the person is indifferent only when the immediate amount equals
  the full later reward (**no discounting**).
- y = 0: any immediate amount, however small, is preferred (**complete
  discounting**).

Real titration data pile up at both ends, and the example data included
with this package (`data(dd_ip)`) reproduce the pattern: a meaningful
fraction of observations sit *exactly* at 0 or 1. That single fact rules
out the two error distributions people reach for first:

- **Gaussian / nonlinear least squares** treats residuals as unbounded
  and homoscedastic. It can predict an IP below 0 or above 1, and it
  assumes the noise is the same size near the bounds as in the middle,
  which is never true for a proportion.
- **Standard beta regression** is the natural model for a proportion,
  but its density is *undefined* at exactly 0 and 1 (the log-likelihood
  is -\infty). Every boundary observation has to be discarded or nudged,
  which biases the fit.

The **Scale-Location-Truncated beta (SLT-beta)** distribution (Kim,
Koffarnus & Franck, 2024; Kim, Kaplan, Koffarnus & Franck, 2025) keeps
the beta’s natural fit to proportion data but assigns *finite
probability* to the endpoints, so 0 and 1 are modeled rather than thrown
away.

## What the SLT-beta density is

Write the mean of the IP at delay D as a discounting function of a
single rate k, using either Mazur’s hyperbola \mu = 1/(1 + kD) or the
exponential \mu = e^{-kD}, and let \phi be a precision parameter (larger
\phi gives a tighter spread around the curve). With shape parameters a =
\mu\phi and b = (1-\mu)\phi, and two *fixed* micro-constants s =
1{+}10^{-7} and l = 10^{-8}, the log-density of an observed IP y \in
\[0,1\] is

\log f(y \mid \mu, \phi) = \underbrace{\log\Gamma(a{+}b) -
\log\Gamma(a) - \log\Gamma(b)}\_{-\log B(a,b)} +
(a{-}1)\log\\\big(\tfrac{y}{s}{+}l\big) + (b{-}1)\log\\\big(1 -
\big(\tfrac{y}{s}{+}l\big)\big) - \log s - \log Z,

where the truncation normalizer is Z =
\mathrm{pbeta}(\tfrac{1}{s}{+}l,\\a,\\b) - \mathrm{pbeta}(l,\\a,\\b).

The intuition is a *microscopic stretch*. The data window y\in\[0,1\] is
mapped just inside the open interval where the ordinary beta density is
finite, so the kernel is never evaluated at exactly 0 or 1. As s\to 1
and l\to 0 the SLT-beta collapses to the ordinary beta.

``` r

slt_logpdf <- function(y, mu, phi, s = 1.0000001, l = 1e-8) {
  a <- mu * phi; b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(a) - lgamma(b) +
    (a - 1) * log(y / s + l) + (b - 1) * log(1 - (y / s + l)) -
    log(s) - log(pbeta(1 / s + l, a, b) - pbeta(l, a, b))
}
yy <- seq(0, 1, length.out = 400)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (p in list(c(mu = 0.5, phi = 4), c(mu = 0.15, phi = 6))) {
  plot(yy, exp(slt_logpdf(yy, p["mu"], p["phi"])), type = "l", lwd = 2,
       xlab = "indifference point", ylab = "density",
       main = sprintf("mu = %.2f, phi = %g", p["mu"], p["phi"]))
  lines(yy, dbeta(yy, p["mu"] * p["phi"], (1 - p["mu"]) * p["phi"]), lty = 2)
  points(c(0, 1), exp(slt_logpdf(c(0, 1), p["mu"], p["phi"])), pch = 16)
}
```

![SLT-beta (solid) vs. ordinary beta (dashed). The SLT-beta carries
finite density to y = 0 and y = 1; the ordinary beta diverges or
vanishes there.](sltb-discounting_files/figure-html/density-1.png)

SLT-beta (solid) vs. ordinary beta (dashed). The SLT-beta carries finite
density to y = 0 and y = 1; the ordinary beta diverges or vanishes
there.

``` r

par(op)
```

Because the variance of a beta is \mu(1-\mu)/(\phi+1), the SLT-beta
gives less noise near the floor and ceiling and more in the middle. That
is the heteroscedastic pattern real IP data show, and it is what a
Gaussian model gets wrong.

## The boundary problem, in the example data

The package includes `dd_ip`, a synthetic delay-discounting dataset (100
subjects across six delays, in the `id`, `x` (delay), `y` (indifference
point) long format) that reproduces the boundary pattern of real
titration data. A handful of synthetic values fall just above 1 and are
clamped to the ceiling.

``` r

data(dd_ip)
dd_ip$id <- factor(dd_ip$id)
dd_ip$y  <- pmin(dd_ip$y, 1)            # clamp a few synthetic values to the [0, 1] ceiling
c(observations        = nrow(dd_ip),
  subjects            = nlevels(dd_ip$id),
  at_floor_or_ceiling = sum(dd_ip$y %in% c(0, 1)))
#>        observations            subjects at_floor_or_ceiling 
#>                 600                 100                 129
```

Many observations sit at exactly 0 and 1. Counting subjects with at
least one boundary value, and asking which error models can even
evaluate them:

``` r

beta_logpdf <- function(y, mu, phi) {          # ordinary beta: -Inf at y in {0,1}
  a <- mu * phi; b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(a) - lgamma(b) + (a - 1) * log(y) + (b - 1) * log(1 - y)
}
ids <- levels(dd_ip$id)
yof <- function(i) dd_ip$y[dd_ip$id == i]
bnd     <- vapply(ids, function(i) any(yof(i) %in% c(0, 1)), logical(1))
beta_ok <- vapply(ids, function(i) all(is.finite(beta_logpdf(yof(i), 0.5, 8))), logical(1))
slt_ok  <- vapply(ids, function(i) all(is.finite(slt_logpdf(yof(i), 0.5, 8))), logical(1))
data.frame(
  subjects                       = length(ids),
  with_boundary_IP               = sum(bnd),
  evaluable_under_ordinary_beta  = sum(beta_ok),
  evaluable_under_SLT_beta       = sum(slt_ok)
)
#>   subjects with_boundary_IP evaluable_under_ordinary_beta
#> 1      100               76                            24
#>   evaluable_under_SLT_beta
#> 1                      100
```

Every subject with a boundary value is lost to ordinary beta regression
and retained by the SLT-beta. Fitting each subject individually
(estimating \log k and \log\phi so that both stay positive) and
comparing the recovered rate to nonlinear least squares confirms that
the SLT-beta agrees with NLS on the bulk of subjects *and* succeeds on
the boundary subjects that NLS-plus-beta could not model:

``` r

fit_slt <- function(y, D) {                      # per-subject SLT MLE (log k, log phi)
  nll <- function(th) {
    mu <- pmin(pmax(1 / (1 + exp(th[1]) * D), 1e-6), 1 - 1e-6)
    v <- -sum(slt_logpdf(y, mu, exp(th[2]))); if (is.finite(v)) v else 1e10
  }
  o <- optim(c(log(0.01), log(8)), nll, control = list(reltol = 1e-11))
  exp(o$par[1])
}
nls_k <- function(y, D) {                        # multi-start NLS (single start is fragile)
  for (st in c(exp(-10), 1e-3, 1e-2, 1e-1)) {
    m <- tryCatch(nls(y ~ 1 / (1 + k * D), start = list(k = st)), error = function(e) NULL)
    if (!is.null(m)) return(coef(m)[["k"]])
  }
  NA_real_
}
xof <- function(i) dd_ip$x[dd_ip$id == i]
k_slt <- vapply(ids, function(i) fit_slt(yof(i), xof(i)), numeric(1))
k_nls <- vapply(ids, function(i) nls_k(yof(i), xof(i)), numeric(1))
sprintf("cor(log k): SLT vs NLS = %.3f   |   median k: SLT = %.4f, NLS = %.4f",
        cor(log(k_slt), log(k_nls), use = "complete.obs"),
        median(k_slt), median(k_nls, na.rm = TRUE))
#> [1] "cor(log k): SLT vs NLS = 0.978   |   median k: SLT = 0.5328, NLS = 0.5038"
```

One caveat, which the next section resolves, is that fitting *one
subject at a time* can let the SLT-beta likelihood run off to a
degenerate solution. For a person whose indifference points are almost
all 0s and 1s, a fit with near-zero precision that places all its mass
on the endpoints can score a higher likelihood than any sensible
discounting curve. The fix is to fit people jointly rather than in
isolation.

## Fitting the mixed-effects model

A mixed-effects model shares information across people. The subject
random intercept on \log k acts as a prior that pulls each person’s
estimate toward the population; the degenerate per-subject solutions
described above are thus regularized away, which is one of the main
reasons to prefer
[`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)
over fitting each subject separately.

The following chunks require TMB (a compiled C++ template;
`has_tmb = TRUE` on this build):

``` r

fit <- fit_dd_tmb(
  data           = dd_ip,            # columns: id, x (delay), y (indifference point)
  equation       = "mazur",          # or "exponential"
  family         = "sltb",           # default; "gaussian" reproduces the NLS-style fit
  random_effects = k ~ 1             # random intercept on log k (fit internally on log scale)
)
#> Fitting TMB mixed-effects discounting model (mazur, sltb)...
#>   Subjects: 100, Observations: 600
#>   Multi-start: best NLL = -2245.78 (start set 1 of 3)
#>   Converged (NLL = -2245.78). Done.
summary(fit)
#> 
#> TMB Mixed-Effects Discounting Model Summary
#> ================================================== 
#> 
#> Call:
#> fit_dd_tmb(data = dd_ip, equation = "mazur", family = "sltb", 
#>     random_effects = k ~ 1) 
#> 
#> Equation: mazur 
#> Family: sltb 
#> Backend: TMB_mixed 
#> Convergence: Yes 
#> Subjects: 100  Observations: 600 
#> 
#> --- Fixed Effects (k) ---
#>           term estimate std.error statistic p.value
#>  k:(Intercept)   0.3655     0.042   -8.7651  <2e-16
#> 
#> --- Variance Components ---
#>                Component Estimate   Scale
#>  sigma_u (log10-k RE SD)   0.4635   log10
#>          phi (precision)   8.9213 natural
#> 
#> --- Fit Statistics ---
#> Log-likelihood: 2245.78 
#> AIC: -4485.56   BIC: -4472.37
tidy(fit)                            # population log k, precision; back-transformed k
#> # A tibble: 3 × 8
#>   term           estimate std.error statistic   p.value component estimate_scale
#>   <chr>             <dbl>     <dbl>     <dbl>     <dbl> <chr>     <chr>         
#> 1 k:(Intercept)     0.366    0.0420     -8.77  1.87e-18 fixed     natural       
#> 2 sigma_u (log1…    0.464   NA          NA    NA        variance  log10         
#> 3 phi (precisio…    8.92    NA          NA    NA        variance  natural       
#> # ℹ 1 more variable: term_display <chr>
ranef(fit)                           # per-subject discount rates
#>       id         u_i           k
#> 1     P1 -0.38236413 0.243029668
#> 2    P10  0.89627486 0.951404759
#> 3   P100  0.11224827 0.412035024
#> 4    P11  0.04970528 0.385427481
#> 5    P12 -0.56628846 0.199711093
#> 6    P13  0.35044270 0.531309819
#> 7    P14 -0.68745303 0.175483796
#> 8    P15  0.37270221 0.544084203
#> 9    P16  0.35772003 0.535452818
#> 10   P17 -0.83495981 0.149920517
#> 11   P18  0.56670341 0.669259120
#> 12   P19 -0.72473673 0.168637582
#> 13    P2  0.14790845 0.428020146
#> 14   P20  0.57456262 0.674896835
#> 15   P21  0.57828154 0.677581086
#> 16   P22 -0.68844168 0.175298717
#> 17   P23 -1.12811078 0.109640892
#> 18   P24  0.59038621 0.686392179
#> 19   P25  0.50477078 0.626449074
#> 20   P26 -2.37264321 0.029045339
#> 21   P27  0.47191776 0.604862885
#> 22   P28  0.77873789 0.839231362
#> 23   P29 -0.13457306 0.316607700
#> 24    P3  0.09144991 0.402988980
#> 25   P30 -0.14907296 0.311745449
#> 26   P31  1.27264370 1.421768897
#> 27   P32  0.65210774 0.733133391
#> 28   P33  0.97853359 1.038713793
#> 29   P34  0.89342397 0.948514125
#> 30   P35  1.13487002 1.227339709
#> 31   P36  0.40341740 0.562216995
#> 32   P37 -0.27346268 0.272985738
#> 33   P38  1.05166888 1.123045514
#> 34   P39 -1.22037768 0.099357979
#> 35    P4  0.03818552 0.380717427
#> 36   P40 -0.09825494 0.329121714
#> 37   P41  0.36177442 0.537774988
#> 38   P42  0.35603175 0.534488810
#> 39   P43 -1.47320898 0.075858492
#> 40   P44 -2.59810504 0.022833099
#> 41   P45  0.93522205 0.991788511
#> 42   P46 -0.58244907 0.196295809
#> 43   P47  0.55310528 0.659615664
#> 44   P48 -1.44965351 0.077789899
#> 45   P49 -1.89016940 0.048609875
#> 46    P5  0.88540814 0.940433541
#> 47   P50  0.51567942 0.633785660
#> 48   P51  0.56939444 0.671184178
#> 49   P52  0.69095671 0.764172191
#> 50   P53  0.34820784 0.530043957
#> 51   P54  0.15132126 0.429582122
#> 52   P55 -2.57732996 0.023345062
#> 53   P56 -0.38855222 0.241429780
#> 54   P57  0.64995458 0.731450457
#> 55   P58 -0.58170758 0.196451225
#> 56   P59  0.51883871 0.635926437
#> 57    P6 -0.17586955 0.302955410
#> 58   P60  0.32237455 0.515628609
#> 59   P61 -0.04535661 0.348238886
#> 60   P62  0.32261646 0.515761763
#> 61   P63  0.17996629 0.442919138
#> 62   P64  0.27369453 0.489521365
#> 63   P65  0.33490743 0.522572503
#> 64   P66 -0.65291374 0.182073827
#> 65   P67 -1.34451942 0.087027706
#> 66   P68 -0.42210011 0.232937763
#> 67   P69  0.30424960 0.505749296
#> 68    P7 -0.18982194 0.298477200
#> 69   P70 -0.47443683 0.220282249
#> 70   P71 -3.62279149 0.007648515
#> 71   P72  1.01860430 1.084102768
#> 72   P73  0.39635484 0.557994813
#> 73   P74 -0.51052301 0.211959045
#> 74   P75  1.10465286 1.188386713
#> 75   P76  0.08854337 0.401740723
#> 76   P77 -0.48818519 0.217073359
#> 77   P78  0.51163505 0.631055656
#> 78   P79 -2.66282840 0.021308982
#> 79    P8  0.14380628 0.426150176
#> 80   P80  0.30080220 0.503891771
#> 81   P81  1.14401118 1.239373246
#> 82   P82  0.31114108 0.509483116
#> 83   P83  0.05790995 0.388817592
#> 84   P84  0.18431222 0.444978451
#> 85   P85 -1.12244959 0.110305401
#> 86   P86  0.92938721 0.985631029
#> 87   P87  0.09506892 0.404548636
#> 88   P88  0.31350294 0.510769108
#> 89   P89 -0.43901133 0.228770909
#> 90    P9  0.41635298 0.570033249
#> 91   P90  0.73021185 0.796870464
#> 92   P91 -0.13013190 0.318112069
#> 93   P92  0.71785045 0.786425654
#> 94   P93  0.40033474 0.560370184
#> 95   P94  1.04563507 1.115836120
#> 96   P95 -1.07720764 0.115762642
#> 97   P96 -1.87521275 0.049392110
#> 98   P97  0.54039046 0.650724374
#> 99   P98  0.73834764 0.803820416
#> 100  P99  0.36389053 0.538990998
```

Group contrasts use the same estimated-marginal-means surface computed
on the `log k` scale and back-transformed, so a between-group comparison
reads as a **ratio of discount rates**. The example below adds a
synthetic two-group factor for illustration:

``` r

set.seed(2)
grp <- data.frame(
  id    = levels(dd_ip$id),
  group = factor(sample(c("A", "B"), nlevels(dd_ip$id), replace = TRUE))
)
dd_grp <- merge(dd_ip, grp, by = "id")
fit_grp <- fit_dd_tmb(dd_grp, equation = "mazur", factors = "group")
#> Fitting TMB mixed-effects discounting model (mazur, sltb)...
#>   Subjects: 100, Observations: 600
#>   Multi-start: best NLL = -2247.59 (start set 3 of 3)
#>   Converged (NLL = -2247.59). Done.
get_dd_param_emms(fit_grp)                                        # EMM of k per group
#> # A tibble: 2 × 6
#>   level       k  k_log std.error conf.low conf.high
#>   <chr>   <dbl>  <dbl>     <dbl>    <dbl>     <dbl>
#> 1 group=A 0.455 -0.788     0.161    0.332     0.623
#> 2 group=B 0.296 -1.22      0.157    0.217     0.403
get_dd_comparisons(fit_grp, compare_specs = ~group, adjust = "holm")  # ratio-of-k contrasts
#> $k
#> $k$emmeans
#> # A tibble: 2 × 6
#>   level       k  k_log std.error conf.low conf.high
#>   <chr>   <dbl>  <dbl>     <dbl>    <dbl>     <dbl>
#> 1 group=A 0.455 -0.788     0.161    0.332     0.623
#> 2 group=B 0.296 -1.22      0.157    0.217     0.403
#> 
#> $k$contrasts_log10
#> # A tibble: 1 × 8
#>   contrast         estimate std.error statistic    df conf.low conf.high p.value
#>   <chr>               <dbl>     <dbl>     <dbl> <dbl>    <dbl>     <dbl>   <dbl>
#> 1 group=A - group…    0.187    0.0974      1.92   Inf -0.00408     0.378  0.0551
#> 
#> $k$contrasts_ratio
#> # A tibble: 1 × 5
#>   contrast          ratio conf.low conf.high p.value
#>   <chr>             <dbl>    <dbl>     <dbl>   <dbl>
#> 1 group=A - group=B  1.54    0.991      2.39  0.0551
#> 
#> 
#> attr(,"class")
#> [1] "beezdiscounting_comparison"
#> attr(,"backend")
#> [1] "tmb"
#> attr(,"adjustment_method")
#> [1] "holm"
#> attr(,"compare_specs_used")
#> [1] "~group"
#> attr(,"contrast_type_used")
#> [1] "pairwise"
#> attr(,"contrast_by_used")
#> [1] "NULL"
```

## Choosing a family

|  | `family = "sltb"` (default) | `family = "gaussian"` |
|----|----|----|
| Respects \[0,1\] | yes (incl. boundary mass at 0 and 1) | no (can predict outside) |
| Variance | shrinks near bounds, grows mid-range | constant |
| Boundary data | modeled | only via the bounded mean |
| Matches legacy [`fit_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd.md) | — | yes (NLS-style) |

Use `"sltb"` for inference on real IP data; use `"gaussian"` when
continuity with a classic least-squares discounting analysis is required
or as a quick baseline.

## How much to trust it

The SLT-beta density was validated before release: it integrates to one,
has finite density at both boundaries, reproduces the published moments,
and recovers known discount rates from simulated data (rank correlation
0.99 for Mazur, 0.95 for the exponential). On the bundled example data,
every boundary subject is fit successfully where ordinary beta
regression cannot be evaluated at all. The full Monte-Carlo study lives
in the package’s development repository on GitHub.

## References

- Kim, M., Koffarnus, M. N., & Franck, C. T. (2024). *Thinking Inside
  the Bounds: Improved Error Distributions for Indifference Point Data
  Analysis and Simulation via Beta Regression Using Common Discounting
  Functions.*
- Kim, M., Kaplan, B. A., Koffarnus, M. N., & Franck, C. T. (2025).
  *Scale-Location-Truncated Beta Regression: Expanding Beta Regression
  to Accommodate 0 and 1.* arXiv:2509.13167.
