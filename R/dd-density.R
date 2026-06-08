# R/dd-density.R

#' SLT-beta and Gaussian log-densities (pure R reference)
#'
#' Pure-R implementations of the two observation-family log-densities used by
#' the mixed-effects discounting model. `.dd_slt_logpdf()` is transcribed
#' verbatim from the verified reference
#' `dev/sltb-verification/verify_sltb.R::slt_logpdf` and is the source of truth
#' against which the compiled TMB template (`src/MixedDiscounting.h`) is
#' cross-checked to 1e-8. See Kim, Koffarnus & Franck (2024) and
#' Kim, Kaplan, Koffarnus & Franck (2025; arXiv:2509.13167).
#'
#' The scale-location-truncation (SLT) constants are the *reference-code*
#' values `s_slt = 1.0000001` (= 1 + 1e-7) and `l = 1e-8`. The truncation
#' normalizer `Z = pbeta(1/s_slt + l, a, b) - pbeta(l, a, b)` is load-bearing
#' (NOT approx 1 at small shapes) and the boundaries are finite at y = 0 and
#' y = 1 by construction.
#'
#' @param y Numeric vector of indifference proportions in `[0, 1]`.
#' @param mu Numeric (scalar or vector recycled against `y`); the mean on the
#'   identity link (a discounting function value).
#' @param phi Numeric (scalar or vector); SLT-beta precision (> 0).
#' @param s_slt,l SLT scale/location constants; defaults are the verified values.
#'   (`s_slt` was previously named `s`; renamed so `s` denotes the discounting
#'   exponent package-wide.)
#' @return Numeric vector of log-density values, the same length as the
#'   recycled inputs.
#' @keywords internal
#' @noRd
.dd_slt_logpdf <- function(y, mu, phi, s_slt = 1.0000001, l = 1e-8) {
  a <- mu * phi
  b <- (1 - mu) * phi
  lgamma(a + b) - lgamma(b) - lgamma(a) +
    (a - 1) * log(y / s_slt + l) + (b - 1) * log(1 - (y / s_slt + l)) -
    log(s_slt) - log(stats::pbeta(1 / s_slt + l, a, b) - stats::pbeta(l, a, b))
}

# .dd_gaussian_logpdf(y, mu, sigma_e)
# Internal Gaussian log-density (not exported; no Rd generated).
# y       - Numeric vector of observed indifference proportions.
# mu      - Numeric (scalar or vector); predicted mean on the identity link.
# sigma_e - Numeric (> 0); Gaussian residual standard deviation.
.dd_gaussian_logpdf <- function(y, mu, sigma_e) {
  stats::dnorm(y, mean = mu, sd = sigma_e, log = TRUE)
}
