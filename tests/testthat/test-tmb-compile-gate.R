# tests/testthat/test-tmb-compile-gate.R
# Compile gate: the C++ MixedDiscounting -nll must equal the verified R SLT
# log-density (and Gaussian density) to 1e-8, isolating the observation
# term by zeroing the random effects.
#
# The R-side density is the canonical .dd_slt_logpdf from R/dd-density.R
# (Task F.1), made available by the `devtools::load_all()` in the run
# commands below. Do NOT inline a local copy here — there is a single source
# of truth for the SLT density. (Gaussian compares against dnorm(log=TRUE).)
#
# *** Marginal vs conditional note ***
# MakeADFun is built WITHOUT `random =` so that obj$fn(par) returns the full
# JOINT (conditional) negative log-likelihood:
#   sum_i  -log f_obs(y_i | mu_i, aux)  +  sum_j  -dnorm(u_j, 0, 1, log=TRUE)
# Using random = "u" would return the Laplace-marginal nll (u integrated out),
# which does NOT equal the R conditional density sum even at u = 0.  The
# RE prior is accounted for explicitly on the R side (re_prior below).

# Build a single-subject, intercept-only TMB data/parameter set and the
# matching R-side prediction for a given equation/family/s. The parameter list
# always includes log_s (the C++ template now declares it); 1-parameter
# equations map it (held at s = 1) via .gate_nll(has_s = FALSE).
.gate_setup <- function(eqn_type, family, beta0, aux, s = 1,
                        delays = c(7, 30, 180, 365, 730, 1460, 2920),
                        y = c(0.95, 0.80, 0.50, 0.35, 0.20, 0.08, 0.00)) {
  n_obs <- length(delays)
  k <- exp(beta0)
  mu_raw <- switch(as.character(eqn_type),
    "0" = 1 / (1 + k * delays),
    "1" = exp(-k * delays),
    "2" = (1 + k * delays)^(-s),
    "3" = ifelse(delays > 0, 1 / (1 + k * delays^s), 1)
  )
  mu <- pmin(pmax(mu_raw, 1e-6), 1 - 1e-6)
  list(
    data = list(
      model = "MixedDiscounting",
      y = y, x = delays,
      subject_id = rep(0L, n_obs),
      X = matrix(1, nrow = n_obs, ncol = 1L),
      eqn_type = eqn_type, family = family,
      n_obs = n_obs, n_subjects = 1L,
      n_re = 1L
    ),
    parameters = list(
      beta_k = beta0,
      log_sigma_u = log(0.5),     # arbitrary; RE drops out because u = 0
      log_aux = log(aux),
      log_s = log(s),
      u = matrix(0, nrow = 1L, ncol = 1L),
      log_sd_re = c(log(0.5), log(0.5)),   # mapped out when n_re == 1
      cor_re = 0,                          # mapped out when n_re == 1
      b = matrix(0, nrow = 1L, ncol = 2L)  # mapped out when n_re == 1
    ),
    mu = mu, y = y
  )
}

# Build the joint (conditional) -nll at the setup's true values. For
# 1-parameter equations log_s is map-fixed; obj$par then excludes it, so
# obj$fn(obj$par) evaluates at exactly the intended point in both cases.
.gate_nll <- function(g, has_s) {
  map <- list(log_sd_re = factor(c(NA, NA)),
              cor_re = factor(NA),
              b = factor(rep(NA, length(g$parameters$b))))
  if (!isTRUE(has_s)) map$log_s <- factor(NA)
  obj <- TMB::MakeADFun(data = g$data, parameters = g$parameters,
                        map = map, DLL = "beezdiscounting", silent = TRUE)
  as.numeric(obj$fn(obj$par))
}

describe("MixedDiscounting compile gate", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("compiles and links the MixedDiscounting template", {
    expect_true("beezdiscounting" %in% names(getLoadedDLLs()))
  })

  it("C++ -nll matches the R SLT log-density to 1e-8 (mazur, sltb)", {
    g <- .gate_setup(eqn_type = 0L, family = 0L, beta0 = log(0.01), aux = 8)
    cpp_nll <- .gate_nll(g, has_s = FALSE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-.dd_slt_logpdf(g$y, g$mu, 8)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches the R SLT log-density to 1e-8 (exponential, sltb)", {
    g <- .gate_setup(eqn_type = 1L, family = 0L, beta0 = log(2e-4), aux = 8)
    cpp_nll <- .gate_nll(g, has_s = FALSE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-.dd_slt_logpdf(g$y, g$mu, 8)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches the R Gaussian density to 1e-8 (mazur, gaussian)", {
    g <- .gate_setup(eqn_type = 0L, family = 1L, beta0 = log(0.01), aux = 0.1)
    cpp_nll <- .gate_nll(g, has_s = FALSE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-dnorm(g$y, g$mu, 0.1, log = TRUE)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches the R Gaussian density to 1e-8 (exponential, gaussian)", {
    g <- .gate_setup(eqn_type = 1L, family = 1L, beta0 = log(2e-4), aux = 0.1)
    cpp_nll <- .gate_nll(g, has_s = FALSE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-dnorm(g$y, g$mu, 0.1, log = TRUE)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches R SLT density at s != 1 (green-myerson, sltb)", {
    g <- .gate_setup(eqn_type = 2L, family = 0L, beta0 = log(0.01), aux = 8,
                     s = 0.7)
    cpp_nll <- .gate_nll(g, has_s = TRUE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-.dd_slt_logpdf(g$y, g$mu, 8)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches R Gaussian density at s != 1 (green-myerson, gaussian)", {
    g <- .gate_setup(eqn_type = 2L, family = 1L, beta0 = log(0.01), aux = 0.1,
                     s = 0.7)
    cpp_nll <- .gate_nll(g, has_s = TRUE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-dnorm(g$y, g$mu, 0.1, log = TRUE)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches R SLT density at s != 1 (rachlin, sltb)", {
    g <- .gate_setup(eqn_type = 3L, family = 0L, beta0 = log(0.01), aux = 8,
                     s = 1.3)
    cpp_nll <- .gate_nll(g, has_s = TRUE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-.dd_slt_logpdf(g$y, g$mu, 8)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches R Gaussian density at s != 1 (rachlin, gaussian)", {
    g <- .gate_setup(eqn_type = 3L, family = 1L, beta0 = log(0.01), aux = 0.1,
                     s = 1.3)
    cpp_nll <- .gate_nll(g, has_s = TRUE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-dnorm(g$y, g$mu, 0.1, log = TRUE)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("rachlin handles an x = 0 row (mu = 1; pow(0,s) AD guard), sltb", {
    g <- .gate_setup(eqn_type = 3L, family = 0L, beta0 = log(0.01), aux = 8,
                     s = 1.3,
                     delays = c(0, 7, 30, 180, 365, 730, 1460),
                     y = c(1.00, 0.95, 0.80, 0.50, 0.35, 0.20, 0.08))
    cpp_nll <- .gate_nll(g, has_s = TRUE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-.dd_slt_logpdf(g$y, g$mu, 8)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
    expect_true(is.finite(cpp_nll))
  })

  it("rachlin handles an x = 0 row (mu = 1; pow(0,s) AD guard), gaussian", {
    g <- .gate_setup(eqn_type = 3L, family = 1L, beta0 = log(0.01), aux = 0.1,
                     s = 1.3,
                     delays = c(0, 7, 30, 180, 365, 730, 1460),
                     y = c(1.00, 0.95, 0.80, 0.50, 0.35, 0.20, 0.08))
    cpp_nll <- .gate_nll(g, has_s = TRUE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-dnorm(g$y, g$mu, 0.1, log = TRUE)) + re_prior
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
    expect_true(is.finite(cpp_nll))
  })
})

describe("MixedDiscounting 2-RE (k, phi) compile gate", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  # R-side reference for the n_re == 2 conditional (joint) -nll. b is the
  # n_subjects x 2 standardized deviate matrix; re = (L b')'.
  .gate2_nll_r <- function(y, x, subject_id, beta0, log_aux, log_sd_re, cor_re,
                           b, eqn_type = 0L, s = 1) {
    sd <- exp(log_sd_re); rho <- tanh(cor_re)
    Sigma <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                      rho * sd[1] * sd[2], sd[2]^2), 2L)
    L <- t(chol(Sigma))
    re <- t(L %*% t(b))                         # n_subjects x 2
    k <- exp(beta0 + re[subject_id + 1L, 1L])
    phi <- pmax(exp(log_aux + re[subject_id + 1L, 2L]), 0.1)  # per-subject floor
    mu_raw <- switch(as.character(eqn_type),
      "0" = 1 / (1 + k * x), "1" = exp(-k * x),
      "2" = (1 + k * x)^(-s), "3" = ifelse(x > 0, 1 / (1 + k * x^s), 1))
    mu <- pmin(pmax(mu_raw, 1e-6), 1 - 1e-6)
    re_prior <- -sum(dnorm(as.numeric(b), 0, 1, log = TRUE))
    sum(-.dd_slt_logpdf(y, mu, phi)) + re_prior
  }

  .gate2_cpp_nll <- function(data, parameters, has_s) {
    map <- list(log_sigma_u = factor(NA),
                u = factor(rep(NA, length(parameters$u))))
    if (!isTRUE(has_s)) map$log_s <- factor(NA)
    obj <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
                          DLL = "beezdiscounting", silent = TRUE)
    as.numeric(obj$fn(obj$par))
  }

  it("C++ 2-RE -nll matches the R reference (mazur, sltb, pdSymm)", {
    delays <- c(7, 30, 180, 365, 730)
    y1 <- c(0.95, 0.80, 0.50, 0.30, 0.12)
    y2 <- c(0.90, 0.70, 0.45, 0.25, 0.10)
    data <- list(model = "MixedDiscounting", y = c(y1, y2),
                 x = rep(delays, 2), subject_id = rep(c(0L, 1L), each = 5L),
                 X = matrix(1, nrow = 10L, ncol = 1L), eqn_type = 0L,
                 family = 0L, n_obs = 10L, n_subjects = 2L, n_re = 2L)
    b <- matrix(c(0.4, -0.3, -0.6, 0.5), nrow = 2L, byrow = TRUE)  # 2 x 2
    parameters <- list(beta_k = log(0.01), log_sigma_u = log(0.5),
                       log_aux = log(8), log_s = 0,
                       u = matrix(0, 2L, 1L),
                       log_sd_re = c(log(0.5), log(0.4)), cor_re = 0.3, b = b)
    cpp <- .gate2_cpp_nll(data, parameters, has_s = FALSE)
    r <- .gate2_nll_r(data$y, data$x, data$subject_id, beta0 = log(0.01),
                      log_aux = log(8), log_sd_re = c(log(0.5), log(0.4)),
                      cor_re = 0.3, b = b, eqn_type = 0L)
    expect_equal(cpp, r, tolerance = 1e-8)
  })

  it("C++ 2-RE -nll matches the R reference (mazur, sltb, pdDiag rho=0)", {
    delays <- c(7, 30, 180, 365, 730)
    data <- list(model = "MixedDiscounting", y = c(0.95, 0.80, 0.50, 0.30, 0.12,
                 0.90, 0.70, 0.45, 0.25, 0.10), x = rep(delays, 2),
                 subject_id = rep(c(0L, 1L), each = 5L),
                 X = matrix(1, nrow = 10L, ncol = 1L), eqn_type = 0L,
                 family = 0L, n_obs = 10L, n_subjects = 2L, n_re = 2L)
    b <- matrix(c(0.4, -0.3, -0.6, 0.5), nrow = 2L, byrow = TRUE)
    # pdDiag: cor_re mapped out at 0 -> rho = tanh(0) = 0.
    parameters <- list(beta_k = log(0.01), log_sigma_u = log(0.5),
                       log_aux = log(8), log_s = 0, u = matrix(0, 2L, 1L),
                       log_sd_re = c(log(0.5), log(0.4)), cor_re = 0, b = b)
    map <- list(log_sigma_u = factor(NA), u = factor(c(NA, NA)),
                cor_re = factor(NA), log_s = factor(NA))
    obj <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
                          DLL = "beezdiscounting", silent = TRUE)
    cpp <- as.numeric(obj$fn(obj$par))
    r <- .gate2_nll_r(data$y, data$x, data$subject_id, beta0 = log(0.01),
                      log_aux = log(8), log_sd_re = c(log(0.5), log(0.4)),
                      cor_re = 0, b = b, eqn_type = 0L)
    expect_equal(cpp, r, tolerance = 1e-8)
  })
})
