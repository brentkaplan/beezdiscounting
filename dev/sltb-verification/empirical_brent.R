## ============================================================================
## Empirical validation of SLT-beta discounting on real data (Data_brent.csv =
## Jarvis et al. 2019: N=126, delays 1wk..8yr). Demonstrates SLT handles the
## boundary IPs (0 and 1) where standard beta regression fails, recovers
## sensible k, and agrees with NLS on non-boundary structure.
## ============================================================================
options(warn = 1)
S <- 1.0000001; L <- 1e-8

slt_logpdf <- function(ip, mu, phi, s = S, l = L) {
  a <- mu * phi; b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(b) - lgamma(a) +
    (a - 1) * log(ip / s + l) + (b - 1) * log(1 - (ip / s + l)) -
    log(s) - log(pbeta(1 / s + l, a, b) - pbeta(l, a, b))
}
## standard (non-truncated) beta log-density — undefined at ip in {0,1}
beta_logpdf <- function(ip, mu, phi) {
  a <- mu * phi; b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(a) - lgamma(b) + (a - 1) * log(ip) + (b - 1) * log(1 - ip)
}

d <- read.csv("/Users/brent/Downloads/Data_brent.csv", check.names = FALSE)
delays <- c(7, 30, 180, 365, 730, 1460, 2920)  # 1wk,1mo,6mo,1yr,2yr,4yr,8yr (days)
## NOTE: d$PID is NOT unique (78 unique / 126 rows). Use row index as the
## subject unit for per-subject fits (each row = one subject-record).
cat(sprintf("PID uniqueness: %d unique / %d rows (%d duplicated)\n\n",
            length(unique(d$PID)), nrow(d), sum(duplicated(d$PID))))
ids <- seq_len(nrow(d))
long <- data.frame(id = rep(ids, each = 7),
                   delay = rep(delays, times = nrow(d)),
                   IP = as.numeric(t(as.matrix(d[, -1]))))

cat(sprintf("subjects = %d, observations = %d, IP range = [%.3f, %.3f]\n",
            nrow(d), nrow(long), min(long$IP), max(long$IP)))
cat(sprintf("boundary IPs: #(IP==0) = %d, #(IP==1) = %d\n",
            sum(long$IP == 0), sum(long$IP == 1)))
nbound <- sum(tapply(long$IP, long$id, function(y) any(y == 0 | y == 1)))
cat(sprintf("subjects with >=1 boundary IP (0 or 1): %d / %d\n\n", nbound, nrow(d)))

## --- the headline demonstration: finiteness on real data ---
beta_nonfinite <- sapply(ids, function(i) {
  y <- long$IP[long$id == i]; any(!is.finite(beta_logpdf(y, 0.5, 8)))
})
slt_finite_all <- all(sapply(ids, function(i) {
  y <- long$IP[long$id == i]; all(is.finite(slt_logpdf(y, 0.5, 8)))
}))
cat(sprintf("standard beta log-density NON-finite (boundary) on %d subjects\n", sum(beta_nonfinite)))
cat(sprintf("SLT  log-density FINITE on all 126 subjects: %s\n\n", slt_finite_all))

## --- per-subject fits: NLS and SLT (log k, log phi; mu guarded off {0,1}) ---
fit_slt <- function(y, D) {
  nll <- function(th) {
    k <- exp(th[1]); phi <- exp(th[2])
    mu <- pmin(pmax(1 / (1 + k * D), 1e-6), 1 - 1e-6)
    v <- -sum(slt_logpdf(y, mu, phi)); if (is.finite(v)) v else 1e10
  }
  o <- tryCatch(optim(c(log(0.01), log(8)), nll, control = list(reltol = 1e-11)),
                error = function(e) NULL)
  if (is.null(o)) c(NA, NA, NA) else c(exp(o$par[1]), exp(o$par[2]), o$convergence)
}
nls_k <- function(y, D) {
  for (st in c(exp(-10), 0.001, 0.01, 0.1)) {
    m <- tryCatch(nls(y ~ 1 / (1 + k * D), start = list(k = st)), error = function(e) NULL)
    if (!is.null(m)) return(coef(m)[["k"]])
  }
  NA_real_
}

res <- data.frame(id = ids, nls_k = NA_real_, slt_k = NA_real_,
                  slt_phi = NA_real_, slt_conv = NA_real_, has_bound = NA_integer_)
for (j in seq_along(ids)) {
  y <- long$IP[long$id == ids[j]]; D <- long$delay[long$id == ids[j]]
  s <- fit_slt(y, D)
  res$nls_k[j] <- nls_k(y, D)
  res$slt_k[j] <- s[1]; res$slt_phi[j] <- s[2]; res$slt_conv[j] <- s[3]
  res$has_bound[j] <- as.integer(any(y == 0 | y == 1))
}

cat(sprintf("SLT fits converged (conv==0): %d / %d\n",
            sum(res$slt_conv == 0, na.rm = TRUE), nrow(res)))
cat(sprintf("  ...including %d / %d subjects that have boundary IPs\n",
            sum(res$slt_conv == 0 & res$has_bound == 1, na.rm = TRUE), sum(res$has_bound == 1)))

both <- is.finite(res$nls_k) & is.finite(res$slt_k) & res$nls_k > 0 & res$slt_k > 0
cat(sprintf("\nSLT-k vs NLS-k agreement: cor(log k) = %.4f over %d subjects\n",
            cor(log(res$nls_k[both]), log(res$slt_k[both])), sum(both)))
cat(sprintf("median k:  NLS = %.5f,  SLT = %.5f\n",
            median(res$nls_k, na.rm = TRUE), median(res$slt_k, na.rm = TRUE)))
cat(sprintf("SLT phi:  median = %.2f, range = [%.2f, %.2f]\n",
            median(res$slt_phi, na.rm = TRUE), min(res$slt_phi, na.rm = TRUE),
            max(res$slt_phi, na.rm = TRUE)))

## sanity: geometric-mean of subject k (the quantity the TMB population fixed
## effect should back-transform to via exp(beta0))
cat(sprintf("\ngeometric-mean SLT k = %.5f (= exp(mean(log k)); target for exp(beta0_pop))\n",
            exp(mean(log(res$slt_k[res$slt_k > 0]), na.rm = TRUE))))
cat(sprintf("median SLT k        = %.5f (target for exp(beta0_pop) under symmetric RE)\n",
            median(res$slt_k, na.rm = TRUE)))

cat("\n==== EMPIRICAL VALIDATION COMPLETE ====\n")
