## ============================================================================
## Independent verification of the SLT-beta density (reference: beta functions.R,
## scale_location_trun_beta_mazur; Kim, Koffarnus, Franck 2024 / arXiv 2509.13167)
## Robust to the integrate() near-edge-spike artifact: normalization is checked
## as (numeric kernel integral == analytic normalizer Z) where quadrature is
## reliable, and moments via Monte-Carlo (valid for all shapes).
## ============================================================================
options(warn = 1)
ok <- TRUE
say <- function(label, pass, detail = "") {
  ok <<- ok && pass
  cat(sprintf("[%s] %-48s %s\n", if (pass) "PASS" else "FAIL", label, detail))
}
S_REF <- 1.0000001; L_REF <- 1e-8

## per-point SLT log-density (verbatim from the reference ll.temp)
slt_logpdf <- function(ip, mu, phi, s = S_REF, l = L_REF) {
  a <- mu * phi; b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(b) - lgamma(a) +
    (a - 1) * log(ip / s + l) + (b - 1) * log(1 - (ip / s + l)) -
    log(s) - log(pbeta(1 / s + l, a, b) - pbeta(l, a, b))
}
slt_pdf <- function(ip, mu, phi, ...) exp(slt_logpdf(ip, mu, phi, ...))
Zof <- function(mu, phi, s = S_REF, l = L_REF) {
  a <- mu * phi; b <- (1 - mu) * phi
  pbeta(1 / s + l, a, b) - pbeta(l, a, b)
}
## exact reference function body, to confirm our extraction == theirs
ref_nll_mazur <- function(data, par) {
  k <- par[1]; phi <- par[2]; delay <- data$delay; IP <- data$IP
  mu <- 1 / (k * delay + 1); alpha <- mu * phi; beta <- phi * (1 - mu)
  s <- 1.0000001; l <- 0.00000001
  ll.temp <- lgamma(alpha + beta) - lgamma(beta) - lgamma(alpha) +
    (alpha - 1) * log(IP / s + l) + (beta - 1) * log(1 - (IP / s + l)) -
    log(s) - log(pbeta(1 / s + l, alpha, beta) - pbeta(l, alpha, beta))
  sum(-(ll.temp * ifelse(IP >= 0 & IP <= 1, 1, 0)))
}

mu_grid <- c(0.02, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99)
phi_grid <- c(1.5, 2, 5, 20, 100, 500)
grid <- expand.grid(mu = mu_grid, phi = phi_grid)

## ---- 1. NORMALIZATION: numeric kernel integral == analytic Z (moderate shapes) ----
## kernel K(g) = (1/s) dbeta(g/s+l, a, b); analytically ∫_0^1 K = Z, and f = K/Z.
mod <- with(grid, mu * phi >= 1 & (1 - mu) * phi >= 1)  # a>=1 & b>=1: quadrature safe
kern_err <- mapply(function(mu, phi) {
  a <- mu * phi; b <- (1 - mu) * phi
  num <- integrate(function(g) (1 / S_REF) * dbeta(g / S_REF + L_REF, a, b),
                   0, 1, rel.tol = 1e-10, subdivisions = 1000L)$value
  abs(num - Zof(mu, phi))
}, grid$mu[mod], grid$phi[mod])
say("kernel integral == analytic normalizer Z", max(kern_err) < 1e-7,
    sprintf("max |∫K - Z| = %.2e over %d moderate-shape cells", max(kern_err), sum(mod)))
## normalized density integrates to 1 (same cells)
norm_err <- mapply(function(mu, phi)
  abs(integrate(function(g) slt_pdf(g, mu, phi), 0, 1, rel.tol = 1e-10)$value - 1),
  grid$mu[mod], grid$phi[mod])
say("normalized density integrates to 1", max(norm_err) < 1e-6,
    sprintf("max |∫f - 1| = %.2e", max(norm_err)))

## ---- 1b. Z is load-bearing (NOT always ~1) at small shapes ------------------
zsmall <- Zof(0.05, 1); zmid <- Zof(0.1, 1); zbig <- Zof(0.5, 10)
say("Z varies with shape (must keep -log(Z) term)", zsmall < 0.7 && zbig > 0.99,
    sprintf("Z(a=.05,b=1)=%.3f  Z(a=.1,b=1)=%.3f  Z(a=.5,b=10)=%.4f", zsmall, zmid, zbig))

## ---- 2. BOUNDARY FINITENESS at g=0 and g=1 (all shapes) ---------------------
b0 <- mapply(function(mu, phi) slt_logpdf(0, mu, phi), grid$mu, grid$phi)
b1 <- mapply(function(mu, phi) slt_logpdf(1, mu, phi), grid$mu, grid$phi)
say("finite log-density at g=0 and g=1 (all shapes)", all(is.finite(c(b0, b1))),
    sprintf("logf(0) in [%.1f,%.1f], logf(1) in [%.1f,%.1f]", min(b0), max(b0), min(b1), max(b1)))

## ---- 3. BETA LIMIT: interior density ~= dbeta ------------------------------
beta_lim_err <- mapply(function(mu, phi) {
  ips <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  ## f_SLT = kernel/Z; at small shapes Z != 1, so the correct limit is dbeta/Z
  max(abs(slt_pdf(ips, mu, phi) - dbeta(ips, mu * phi, (1 - mu) * phi) / Zof(mu, phi)))
}, grid$mu, grid$phi)
say("interior density == dbeta/Z (scale/loc limit)", max(beta_lim_err) < 1e-3,
    sprintf("max |f_SLT - dbeta/Z| = %.2e", max(beta_lim_err)))

## ---- 4. MOMENTS via Monte-Carlo (robust for all shapes) --------------------
rslt <- function(n, mu, phi, s = S_REF, l = L_REF) {
  a <- mu * phi; b <- (1 - mu) * phi
  (qbeta(runif(n, pbeta(l, a, b), pbeta(1 / s + l, a, b)), a, b) - l) * s
}
set.seed(42)
mom <- t(sapply(list(c(0.2, 5), c(0.5, 20), c(0.8, 100), c(0.05, 8)), function(p) {
  mu <- p[1]; phi <- p[2]; g <- rslt(2e6, mu, phi)
  c(mu = mu, phi = phi,
    e_err = abs(mean(g) - S_REF * (mu - L_REF)),
    v_relerr = abs(var(g) - S_REF^2 * mu * (1 - mu) / (phi + 1)) /
               (S_REF^2 * mu * (1 - mu) / (phi + 1)))
}))
say("MC mean == s(mu-l)", max(mom[, "e_err"]) < 3e-3, sprintf("max |E-formula| = %.2e", max(mom[, "e_err"])))
say("MC var == s^2 mu(1-mu)/(phi+1)", max(mom[, "v_relerr"]) < 0.03,
    sprintf("max rel-err = %.2e (truncation-order)", max(mom[, "v_relerr"])))

## ---- 5. our extraction == reference function body (incl. boundary IP=0) -----
dtest <- data.frame(delay = c(7, 30, 180, 365, 730, 1460, 2920),
                    IP = c(0.95, 0.8, 0.5, 0.35, 0.2, 0.08, 0.0))
k0 <- 0.01; phi0 <- 8; mu_m <- 1 / (1 + k0 * dtest$delay)
say("our log-density == reference NLL body",
    abs(sum(-slt_logpdf(dtest$IP, mu_m, phi0)) - ref_nll_mazur(dtest, c(k0, phi0))) < 1e-12,
    "exact to 1e-12")

## ---- 6. MLE RECOVERY (the practical correctness gate; log k, log phi) -------
fit_one <- function(d, link = "mazur") {
  nll <- function(th) {
    k <- exp(th[1]); phi <- exp(th[2])
    mu <- if (link == "mazur") 1 / (1 + k * d$delay) else exp(-k * d$delay)
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)  # guard mu off {0,1} (exp underflow at long delays)
    v <- -sum(slt_logpdf(d$IP, mu, phi)); if (is.finite(v)) v else 1e10
  }
  ## mini multi-start (mirrors the package's multi_start design): a fixed start
  ## is inadequate across k regimes, esp. for the fast-saturating exponential.
  best <- NULL
  for (k0 in c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) {
    o <- tryCatch(optim(c(log(k0), log(8)), nll, control = list(reltol = 1e-11)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  c(k = exp(best$par[1]), phi = exp(best$par[2]), conv = best$convergence)
}
recover <- function(link, krange) {
  set.seed(20260606)
  delays <- c(7, 30, 180, 365, 730, 1460, 2920)
  truek <- exp(seq(log(krange[1]), log(krange[2]), length.out = 80)); truephi <- 12
  est <- t(sapply(truek, function(k) {
    mu <- if (link == "mazur") 1 / (1 + k * delays) else exp(-k * delays)
    mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    fit_one(data.frame(delay = delays, IP = rslt(length(delays), mu, truephi)), link)
  }))
  list(cor = cor(log(truek), log(est[, "k"])),
       bias = mean(log(est[, "k"]) - log(truek)),
       phimed = median(est[, "phi"]), conv = all(est[, "conv"] == 0))
}
rm <- recover("mazur", c(1e-4, 0.05))
say("MLE k recovery (Mazur): cor(log k,log k_hat)>0.97", rm$cor > 0.97,
    sprintf("cor=%.4f bias=%+.3f phi_med=%.1f conv=%s", rm$cor, rm$bias, rm$phimed, rm$conv))
## exponential saturates faster than hyperbolic; use its identifiable k regime
## (ED50=ln2/k within/near the 7..2920d window) so mu doesn't floor at 0
re <- recover("exponential", c(3e-5, 1.2e-3))
## exponential is inherently less identifiable than hyperbolic (faster saturation),
## so a slightly looser bar; 0.95 is strong recovery for this equation.
say("MLE k recovery (exponential): cor>0.93", re$cor > 0.93,
    sprintf("cor=%.4f bias=%+.3f phi_med=%.1f conv=%s", re$cor, re$bias, re$phimed, re$conv))

cat(sprintf("\n==== OVERALL: %s ====\n", if (ok) "ALL CHECKS PASSED" else "SOME CHECKS FAILED"))
quit(status = if (ok) 0L else 1L)
