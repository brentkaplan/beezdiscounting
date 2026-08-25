.lin_ftest_sim <- function(seed = 11) {
  set.seed(seed)
  lv <- c("EFT", "HIT", "NCC", "WL")
  n <- c(25, 30, 35, 20); T <- 6; mu <- c(-6.5, -6.0, -5.5, -5.8); sigma2 <- 2; g <- 10
  cond <- rep(lv, n)
  theta <- stats::rnorm(sum(n), mu[match(cond, lv)], sqrt(g * sigma2 / T))
  data.frame(
    unit = factor(rep(seq_len(sum(n)), each = T)), condition = factor(rep(cond, each = T), levels = lv),
    y_lin = stats::rnorm(sum(n) * T, rep(theta, each = T), sqrt(sigma2)), log_jac = 1
  )
}
# ORACLE built ONLY from the raw simulated frame `d` (never from the fit object):
# unit means + condition labels aggregated independently, then nested lm ANOVA.
.oracle_F <- function(d, hypothesis) {
  logk <- tapply(d$y_lin, d$unit, mean)
  full <- factor(tapply(as.character(d$condition), d$unit, `[`, 1), levels = levels(d$condition))
  map <- stats::setNames(levels(full), levels(full))
  if (is.null(hypothesis)) hypothesis <- list(levels(full))
  for (set in hypothesis) if (length(set) > 1L) map[set] <- paste(set, collapse = "=")
  dd <- data.frame(logk = as.numeric(logk), full = full, red = factor(map[as.character(full)]))
  # a single-level `red` (all levels merged) is intercept-only; `~ red` errors on contrasts
  fml_red <- if (nlevels(dd$red) > 1L) logk ~ red else logk ~ 1
  a <- anova(lm(fml_red, dd), lm(logk ~ full, dd))
  list(F = a$F[2], df1 = a$Df[2], df2 = a$Res.Df[2], p = a$`Pr(>F)`[2])
}

describe(".dd_lin_ftest()", {
  d <- .lin_ftest_sim()
  re <- .dd_lin_re_mle(d$y_lin, d$unit, d$condition, d$log_jac)

  it("merges levels per hypothesis; NULL collapses everything", {
    m <- .dd_lin_merge_levels(re$unit_condition, list(c("EFT", "NCC")))
    expect_equal(levels(m), c("EFT=NCC", "HIT", "WL"))
    expect_equal(nlevels(.dd_lin_merge_levels(re$unit_condition, NULL)), 1L)
    expect_error(.dd_lin_merge_levels(re$unit_condition, list("nope")), "unknown")
  })
  it("ORACLE all-equal H0: F, df, p equal nested one-way ANOVA on unit means", {
    r <- .dd_lin_ftest(re)
    o <- .oracle_F(d, NULL)
    expect_equal(r$F, o$F, tolerance = 1e-10); expect_equal(r$df1, 3); expect_equal(r$df2, 106)
    expect_equal(r$p_value, o$p, tolerance = 1e-10)
    expect_true(is.na(r$cohens_d))
  })
  it("ORACLE pairwise H0 mu_EFT = mu_NCC (paper Table 3 pattern)", {
    h <- list(c("EFT", "NCC"))
    r <- .dd_lin_ftest(re, h); o <- .oracle_F(d, h)
    expect_equal(r$F, o$F, tolerance = 1e-10); expect_equal(r$df1, 1); expect_equal(r$df2, 106)
    expect_equal(r$p_value, o$p, tolerance = 1e-10)
    expect_equal(r$hypothesis, "EFT = NCC")
  })
  it("ORACLE two simultaneous disjoint constraints: EFT=HIT and NCC=WL (df1 = 2)", {
    h <- list(c("EFT", "HIT"), c("NCC", "WL"))
    r <- .dd_lin_ftest(re, h); o <- .oracle_F(d, h)
    expect_equal(r$F, o$F, tolerance = 1e-10); expect_equal(r$df1, 2); expect_equal(o$df1, 2)
    expect_equal(r$p_value, o$p, tolerance = 1e-10)
  })
  it("ORACLE three-way subset constraint EFT=HIT=NCC leaves WL free (df1 = 2)", {
    h <- list(c("EFT", "HIT", "NCC"))
    r <- .dd_lin_ftest(re, h); o <- .oracle_F(d, h)
    expect_equal(r$F, o$F, tolerance = 1e-10); expect_equal(r$df1, 2)
  })
  it("denominator df is sum(n_c - 1) and uses SSR_{Z|X_full}, not SSE_Z", {
    r <- .dd_lin_ftest(re)
    expect_equal(r$df2, sum(c(25, 30, 35, 20) - 1))
  })
  it("a bare character vector hypothesis behaves exactly like the one-element list", {
    expect_identical(.dd_lin_ftest(re, c("EFT", "NCC")), .dd_lin_ftest(re, list(c("EFT", "NCC"))))
    bare <- .dd_lin_ftest(re, c("EFT", "NCC"))
    expect_equal(bare$hypothesis, "EFT = NCC")
    expect_false(is.na(bare$cohens_d))
  })
  it("Cohen's d for a pair follows Sec 4.3: sqrt(T)(mu1-mu2)/(sigma sqrt(g+1))", {
    r <- .dd_lin_ftest(re, list(c("EFT", "NCC")))
    expect_equal(r$cohens_d, sqrt(6) * (re$mu[["EFT"]] - re$mu[["NCC"]]) / (sqrt(re$sigma2) * sqrt(re$g + 1)))
  })
})
